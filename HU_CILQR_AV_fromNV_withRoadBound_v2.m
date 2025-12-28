%% HU_CILQR_AV_fromNV_withRoadBound_v2_CN.m
% 关键修正：
% 1) 时间尺度对齐：将 nuScenes(通常0.5s) 的 NV 未来轨迹插值到规划步长(0.1s)
% 2) NV 前向轨迹从“当前时刻”开始（而不是从 past 的开头开始）
% 3) 参考车道线选择：距离 + 航向一致性（避免选错/反向车道）
% 4) Stage2 使用发散型 barrier（-log）让安全距离更接近硬约束行为
% 5) 3D 轨迹 z 轴使用真实时间（秒）
% 6) v/a/omega 直接读状态/控制，不用差分

clear; clc; close all;

%% ==================== 0) 输入参数 ====================
jsonPath  = 'ex_list_sample_updated.json';
sampleIdx = 1;

road_bound_showFlag = true;   % 是否画道路边界（用于检查）
doPlotScene = true;           % 是否先画场景做 sanity check

% 数据时间步长（nuScenes labels/past 常见约为 0.5s）
dt_data = 0.5;   % 重要：若你的导出间隔不是0.5s，请改这里
T_s     = 0.1;   % 规划器步长

% 预测时域与仿真时长
N      = 40;      % 0.1s * 40 = 4s
T_sim  = 8.0;     % 期望仿真 8s（会根据 future 长度自动截断）

N_iter  = 8;      % 每个阶段的迭代次数（简化求解器）
N_stage = 2;      % 两阶段：1软约束，2更硬的约束
lambda_damp = 1e-3; % % 预留：以后换回完整 iLQR 时可用

% 车辆长度（双圆模型）
L_AV = 4.5;
L_NV = 4.5;

% 代价权重（来自你 markdown / 论文）
w1 = 2.0; w2 = 0.1; w3 = 1.0; w4 = 3.0;

% 控制输入约束
a_min = -5; a_max = 3;
delta_lim = 0.5;
delta_min = -delta_lim; delta_max = delta_lim;

% 安全距离（此版本不做不确定性膨胀）
s_safe = 2.4;
margin_safe = 1.0;
s_safe_eff = s_safe + margin_safe;

% 速度非负 barrier
v_min = 0.0;
eps_v = 1e-3;

% warm start
U_warm = zeros(2, N);

%% ==================== 1) 从 JSON 读取场景 ====================
scene = get_scene_from_json_withRoadBound(jsonPath, sampleIdx, road_bound_showFlag, doPlotScene);

lanes   = scene.lane_centerlines;     % 车道中心线集合（cell）
road_bd = scene.road_boundaries;      % 道路边界集合（cell）
nv_past = scene.nv_past_xy;           % NV 历史轨迹（这里 NV=ego）
nv_fut  = scene.nv_future_xy;         % NV 未来 GT labels（可为空）
sample  = scene.sample;               % 原始 sample 结构

% NV 当前点（local frame）
p_nv0 = nv_past(end,:);

% NV 航向：用历史最后两点估计
if size(nv_past,1) >= 2
    dp = nv_past(end,:) - nv_past(end-1,:);
    theta0 = atan2(dp(2), dp(1));
else
    theta0 = 0;
end

% NV 速度：优先读 past_traj 第3列，否则用差分估计
v0 = estimate_speed_from_past(sample, nv_past, dt_data);

%% ==================== 2) 构造 NV 的“前向轨迹”（修正：从当前时刻开始） ====================
% 正确的 forward：NV_current + NV_future_labels
nv_forward_raw = p_nv0;
if ~isempty(nv_fut)
    nv_forward_raw = [nv_forward_raw; nv_fut];
end

% 根据 future 的长度决定可仿真时长
t_raw_end = (size(nv_forward_raw,1)-1) * dt_data;
T_sim = min(T_sim, t_raw_end);
K_exec = floor(T_sim / T_s);                 % MPC 执行步数
t_plan = (0:K_exec) * T_s;                   % 规划时刻序列（秒）

% 将 NV 轨迹插值到规划步长
t_raw = (0:size(nv_forward_raw,1)-1) * dt_data;
nv_plan_xy = interp1(t_raw, nv_forward_raw, t_plan, 'linear', 'extrap'); % [(K_exec+1) x 2]

% 若 future 为空，则 NV 位置保持不动
if isempty(nv_fut)
    nv_plan_xy = repmat(p_nv0, numel(t_plan), 1);
end

%% ==================== 3) 构造 AV 初始状态（相对 NV 偏移） ====================
offset_xy = [-2.0, 0];  % 横向偏移，可自行调
p_av0 = p_nv0 + offset_xy;

% 状态 [x;y;v;theta]
S0 = [p_av0(1); p_av0(2); v0; theta0];

%% ==================== 4) 选参考车道线（加入航向一致性） ====================
lane_ref = pick_nearest_lane_with_heading(lanes, p_nv0, theta0);

% 对参考车道线按弧长重采样，方便按“点序列”滚动截取
ds = 0.5;
wpt = resample_polyline_by_arclen(lane_ref, ds);

%% ==================== 5) 道路边界构造 polyshape（用于边界 barrier） ====================
roadPoly = build_road_polyshape_general(road_bd);

%% ==================== 6) MPC 循环 ====================
S_exec = zeros(4, K_exec+1);
U_exec = zeros(2, K_exec);
S_exec(:,1) = S0;
WPT_SEG_LOG = cell(K_exec,1);   % 存每一步的 wpt_seg（用于可视化）


for k0 = 1:K_exec
    S_init = S_exec(:,k0);

    % 取 NV 在 horizon 内的轨迹段（已对齐到 0.1s）
    nv_seg = build_nv_segment_from_plan(nv_plan_xy, k0, N, T_s);

    % 取参考轨迹段（从离当前 AV 最近点开始取 N+1 个点）
    wpt_seg = slice_wpt_for_horizon(wpt, [S_init(1), S_init(2)], N);
    x_wpt = wpt_seg(:,1);  y_wpt = wpt_seg(:,2);
    WPT_SEG_LOG{k0} = wpt_seg;   % 记录本次 MPC 的参考段


    % 求解（简化版：数值梯度下降，占位符结构，先让效果正确）
    [~, U_opt] = solve_HU_CILQR_noSigma_withRoadBound( ...
        S_init, U_warm, x_wpt, y_wpt, v0, ...
        nv_seg, L_AV, L_NV, s_safe_eff, ...
        roadPoly, ...
        w1,w2,w3,w4, a_min,a_max, delta_min,delta_max, v_min, eps_v, ...
        T_s, N, N_stage, N_iter);

    % 执行第一个控制输入
    U0 = U_opt(:,1);
    U_exec(:,k0) = U0;
    S_exec(:,k0+1) = f_dynamics(S_exec(:,k0), U0, T_s);

    % warm start：右移控制序列
    U_warm = [U_opt(:,2:end), U_opt(:,end)];
end




%% ==================== 7) XY 场景绘图 ====================
figure('Name','AV vs NV + lanes + road boundary'); hold on; grid on; axis equal;

% 车道中心线
for l = 1:numel(lanes)
    xy = lanes{l};
    if isempty(xy), continue; end
    plot(xy(:,1), xy(:,2), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
end

% 道路边界（优先画 polyshape 的边界）
if ~isempty(roadPoly)
    [bx,by] = boundary(roadPoly);
    plot(bx, by, '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.8);
end


% ---------- 绘制每一步的 wpt_seg（黑色线段） ----------
% 只画一部分步数，避免太密（可调：比如每5步画一次）
step_draw = 5;

h_wpt = gobjects(0);
for k = 1:step_draw:K_exec
    seg = WPT_SEG_LOG{k};
    if isempty(seg) || size(seg,2)~=2
        continue;
    end

    if isempty(h_wpt)
        h_wpt = plot(seg(:,1), seg(:,2), 'k-', 'LineWidth', 2.0, ...
            'DisplayName','wpt\_seg (reference horizon)');
    else
        plot(seg(:,1), seg(:,2), 'k-', 'LineWidth', 2.0, ...
            'HandleVisibility','off');
    end

    % 可选：在每段起点打点，方便看“参考段起始点是否乱跳”
    plot(seg(1,1), seg(1,2), 'ko', 'MarkerSize', 4, 'HandleVisibility','off');
end



% NV（插值后的轨迹）
plot(nv_plan_xy(:,1), nv_plan_xy(:,2), 'r--', 'LineWidth', 2.0);

% AV 执行轨迹
plot(S_exec(1,:), S_exec(2,:), 'b-', 'LineWidth', 2.5);

xlabel('x_{rel} (m)'); ylabel('y_{rel} (m)');
title(sprintf('Sample #%d: AV (planned) vs NV (Replay)', sampleIdx));
legend({'Lane centerlines','Road boundary','NV (replay)','AV (planned)'}, 'Location','bestoutside');

%% ==================== 8) 3D 轨迹（z=时间秒） ====================
t_exec = (0:K_exec) * T_s;

figure('Name', '3D Trajectories (z=time in seconds)', 'Color', 'w');
hold on; grid on;
plot3(nv_plan_xy(:,1), nv_plan_xy(:,2), t_exec, 'r--', 'LineWidth', 2);
plot3(S_exec(1,:),     S_exec(2,:),     t_exec, 'b-',  'LineWidth', 2.5);
xlabel('x_{rel} (m)'); ylabel('y_{rel} (m)'); zlabel('Time (s)');
title('3D Trajectories: AV vs NV (z = time)');
legend({'NV (replay)','AV (planned)'}, 'Location', 'bestoutside');
view(45,25);
hold off;

%% ==================== 9) 评估曲线（subplot） ====================
% 速度/加速度/角速度：直接来自状态/控制，更稳定、更符合物理意义
v_av = S_exec(3,:);
a_av = U_exec(1,:);
w_av = U_exec(2,:);

% 与 NV 距离
dist_to_nv = sqrt((S_exec(1,:) - nv_plan_xy(:,1)').^2 + (S_exec(2,:) - nv_plan_xy(:,2)').^2);

figure('Name','v / omega / a / distance', 'Color','w');

subplot(4,1,1); grid on; hold on;
plot(t_exec, v_av, 'b-', 'LineWidth', 2);
ylabel('v (m/s)'); title('AV velocity');

subplot(4,1,2); grid on; hold on;
plot(t_exec(1:end-1), w_av, 'r-', 'LineWidth', 2);
ylabel('\omega (rad/s)'); title('AV yaw rate');

subplot(4,1,3); grid on; hold on;
plot(t_exec(1:end-1), a_av, 'g-', 'LineWidth', 2);
ylabel('a (m/s^2)'); title('AV acceleration');

subplot(4,1,4); grid on; hold on;
plot(t_exec, dist_to_nv, 'k-', 'LineWidth', 2);
yline(s_safe_eff,'--','LineWidth',1.5);
xlabel('Time (s)'); ylabel('dist (m)');
title('Distance to NV (and safety threshold)');

%% ========================= 辅助函数 =========================

function out = normalizePolylineList(x)
% 将 jsondecode 后的 “多条 polyline” 统一成 cell{K}，每个元素是 [N x 2] double
    out = {};
    if isempty(x), return; end

    if isnumeric(x) && size(x,2)==2
        out = {double(x)}; return;
    end

    if iscell(x)
        tmp = {};
        for i = 1:numel(x)
            yi = x{i};
            if isempty(yi), continue; end
            if isnumeric(yi) && size(yi,2)==2
                tmp{end+1,1} = double(yi); %#ok<AGROW>
            elseif iscell(yi)
                sub = normalizePolylineList(yi);
                for k = 1:numel(sub)
                    tmp{end+1,1} = sub{k}; %#ok<AGROW>
                end
            end
        end
        out = tmp; return;
    end

    warning('road_boundaries 结构不是 numeric/cell，当前脚本未覆盖该格式。');
end

function v0 = estimate_speed_from_past(m, nv_past_xy, dt_data)
% 优先从 past_traj 第3列 v 取；否则差分估计
    v0 = 5.0;
    try
        if isfield(m,'past_traj') && ~isempty(m.past_traj) && size(m.past_traj,2) >= 3
            v0 = m.past_traj(end,3);
        else
            if size(nv_past_xy,1) >= 2
                dp = nv_past_xy(end,:) - nv_past_xy(end-1,:);
                v0 = norm(dp) / dt_data;
            end
        end
    catch
    end
end

function lane_ref = pick_nearest_lane_with_heading(lanes, p0, theta_av)
% 评分：score = 距离 + gamma*(1 - cos(航向差))
    bestScore = inf; best = [];
    gamma = 8.0; % 航向一致性权重（可调 3~15）

    for i = 1:numel(lanes)
        xy = lanes{i};
        if isempty(xy) || size(xy,1)<3, continue; end

        d = vecnorm(xy - p0, 2, 2);
        [dmin, idx] = min(d);

        idx2 = min(idx+1, size(xy,1));
        dp = xy(idx2,:) - xy(idx,:);
        theta_lane = atan2(dp(2), dp(1));

        score = dmin + gamma*(1 - cos(theta_lane - theta_av));
        if score < bestScore
            bestScore = score;
            best = xy;
        end
    end
    lane_ref = best;
end

function P = resample_polyline_by_arclen(xy, ds)
% 按弧长重采样 polyline
    if isempty(xy) || size(xy,1) < 2
        P = xy; return;
    end
    dseg = vecnorm(diff(xy,1,1),2,2);
    s = [0; cumsum(dseg)];
    s_new = (0:ds:s(end))';
    x_new = interp1(s, xy(:,1), s_new, 'linear', 'extrap');
    y_new = interp1(s, xy(:,2), s_new, 'linear', 'extrap');
    P = [x_new, y_new];
end

function wpt_seg = slice_wpt_for_horizon(wpt, p_now, N)
% 从 wpt 中找最近点作为起始，然后取 N+1 个点（不足就重复末尾）
    d = vecnorm(wpt - p_now,2,2);
    [~,i0] = min(d);
    i1 = i0 + N;

    if i1 <= size(wpt,1)
        wpt_seg = wpt(i0:i1,:);
    else
        tail = repmat(wpt(end,:), i1 - size(wpt,1), 1);
        wpt_seg = [wpt(i0:end,:); tail];
    end
end

function nv_seg = build_nv_segment_from_plan(nv_plan_xy, k0, N, T_s)
% nv_plan_xy: 已对齐到 T_s 的轨迹
    idx0 = k0;
    idx1 = min(k0+N, size(nv_plan_xy,1));
    xy = nv_plan_xy(idx0:idx1,:);

    if size(xy,1) < N+1
        xy = [xy; repmat(xy(end,:), N+1-size(xy,1), 1)];
    end

    v = zeros(N,1);
    th = zeros(N,1);
    for k = 1:N
        dp = xy(k+1,:) - xy(k,:);
        th(k) = atan2(dp(2), dp(1));
        v(k)  = norm(dp) / T_s;
    end

    nv_seg = struct('xy',xy,'v',v,'th',th);
end

function S_next = f_dynamics(S, U, T_s)
% 简化自行车模型（同你原始形式）
    S_next = S + [S(3)*cos(S(4));
                  S(3)*sin(S(4));
                  U(1);
                  U(2)] * T_s;
end

function [S_traj, U_best] = solve_HU_CILQR_noSigma_withRoadBound( ...
    S0, U_init, x_wpt, y_wpt, v_d, ...
    nv_seg, L_AV, L_NV, s_safe_eff, ...
    roadPoly, ...
    w1,w2,w3,w4, a_min,a_max, delta_min,delta_max, v_min, eps_v, ...
    T_s, N, N_stage, N_iter)

% 简化求解器：用数值梯度下降做“可跑占位”
% Stage1：soft（exp barrier）
% Stage2：hard-like（-log barrier）
    U = U_init;
    S_traj = zeros(4,N+1);

    for stage = 1:N_stage
        for it = 1:N_iter %#ok<NASGU>
            % 前向 rollout
            S_traj(:,1) = S0;
            for k = 1:N
                S_traj(:,k+1) = f_dynamics(S_traj(:,k), U(:,k), T_s);
            end

            % 梯度下降更新
            alpha = (stage==1)*0.25 + (stage==2)*0.15;

            for k = 1:N
                u0 = U(:,k);
                g = numeric_grad_u(@(uu) running_cost_one(S_traj(:,k), uu, ...
                        [x_wpt(k);y_wpt(k)], v_d, ...
                        nv_seg.xy(k,:), nv_seg.th(min(k,end)), ...
                        L_AV, L_NV, s_safe_eff, ...
                        roadPoly, stage, v_min, eps_v, ...
                        w1,w2,w3,w4), u0);

                U(:,k) = U(:,k) - alpha*g;

                % 输入限幅
                U(1,k) = min(max(U(1,k), a_min), a_max);
                U(2,k) = min(max(U(2,k), delta_min), delta_max);
            end
        end
    end

    U_best = U;
end

function J = running_cost_one(S, U, wpt, v_d, nv_xy, nv_th, L_AV, L_NV, s_safe_eff, roadPoly, stage, v_min, eps_v, w1,w2,w3,w4)
% running cost：轨迹跟踪 + 控制能量 + 安全 barrier + 道路边界 barrier
    ex = S(1)-wpt(1);
    ey = S(2)-wpt(2);
    ev = S(3)-v_d;
    J = w1*(ex^2 + ey^2) + w2*(ev^2) + w3*(U(1)^2) + w4*(U(2)^2);

    % 速度非负 barrier（Stage2更硬）
    if stage==1
        wv = 10;
        den = (S(3)-v_min+eps_v);
        if den < 0.5
            J = J + wv/den;
        end
    else
        wv = 60;
        den = (S(3)-v_min+eps_v);
        if den <= 0
            J = J + 1e6;
        elseif den < 0.6
            J = J - wv*log(den);
        end
    end

    % 双圆模型最小距离
    dmin = doublecircle_dmin_xytheta(S, nv_xy, nv_th, L_AV, L_NV);
    gap = dmin - s_safe_eff;

    % 安全距离 barrier：Stage1 soft, Stage2 hard-like
    if stage==1
        wobs=35; range=6;
        if gap < range
            J = J + wobs*exp(-gap);
        end
    else
        wobs=120;
        if gap <= 0
            J = J + 1e6;
        elseif gap < 4.0
            J = J - wobs*log(gap);
        end
    end

    % 道路边界 barrier
    if ~isempty(roadPoly)
        [in,on] = isinterior(roadPoly, S(1), S(2));
        d2bd = min_distance_to_polyshape_boundary(roadPoly, [S(1),S(2)]);

        if ~(in || on)
            J = J + (stage==1)*500 + (stage==2)*5000;
            J = J + 2000*(1 + d2bd);
        else
            if stage==1
                if d2bd < 2.0
                    J = J + 10*exp(-(d2bd));
                end
            else
                if d2bd <= 0
                    J = J + 1e6;
                elseif d2bd < 1.5
                    J = J - 40*log(d2bd);
                end
            end
        end
    end
end

function g = numeric_grad_u(fun, u)
% 数值梯度（2维输入）
    eps = 1e-3;
    g = zeros(2,1);
    for i=1:2
        du = zeros(2,1); du(i)=eps;
        g(i) = (fun(u+du)-fun(u-du))/(2*eps);
    end
end

function dmin = doublecircle_dmin_xytheta(S_av, xy_nv, th_nv, L_AV, L_NV)
% 双圆模型：计算 AV 与 NV 的四对圆心距离的最小值
    x = S_av(1); y=S_av(2); th=S_av(4);
    Fe = [x;y] + (L_AV/2)*[cos(th); sin(th)];
    Be = [x;y] - (L_AV/2)*[cos(th); sin(th)];

    Fo = xy_nv(:) + (L_NV/2)*[cos(th_nv); sin(th_nv)];
    Bo = xy_nv(:) - (L_NV/2)*[cos(th_nv); sin(th_nv)];

    dmin = min([norm(Fe-Fo), norm(Fe-Bo), norm(Be-Fo), norm(Be-Bo)]);
end

function roadPoly = build_road_polyshape_general(road_bd)
% 从 road_bd 中挑两条最长边界，拼成一个闭合 polyshape（近似道路区域）
    roadPoly = [];
    try
        if isempty(road_bd), return; end
        if ~iscell(road_bd), return; end

        lens = zeros(numel(road_bd),1);
        for i=1:numel(road_bd)
            xy = road_bd{i};
            if isempty(xy) || size(xy,2)~=2, continue; end
            lens(i) = sum(vecnorm(diff(xy,1,1),2,2));
        end

        [~,ord] = sort(lens,'descend');
        if numel(ord) < 2, return; end

        L = road_bd{ord(1)};
        R = road_bd{ord(2)};

        if size(L,2)==2 && size(R,2)==2 && size(L,1)>2 && size(R,1)>2
            P = [L; flipud(R)];
            roadPoly = polyshape(P(:,1), P(:,2));
        end
    catch
        roadPoly = [];
    end
end

function d = min_distance_to_polyshape_boundary(poly, p)
% 点到 polyshape 边界最小距离（近似：边界点采样）
    d = 0;
    try
        [bx,by] = boundary(poly);
        if isempty(bx), d = 0; return; end
        B = [bx(:), by(:)];
        d = min(vecnorm(B - p,2,2));
    catch
        d = 0;
    end
end
