%% HU_CILQR_AV_fromNV_withRoadBound.m
% - NV: 直接用 JSON 里的 ego 轨迹（past + future GT）
% - AV: 虚构，初始相对 NV 偏移一点点，速度/航向与 NV 一致
% - Ref: 从 lane centerlines 中选“离 NV 当前最近”的那条，作为 x^{wpt}/y^{wpt}
% - Road boundary: 使用 road_boundaries 做 barrier（不改 matrix，只用 mapping 新字段）
% - 不使用 Σ_t / 不膨胀安全距离（按你要求删掉）

clear; clc; close all;

%% ========== 1) 读取场景（函数返回 lanes + road_bd + NV轨迹） ==========
jsonPath = 'ex_list_sample_updated.json';
sampleIdx = 1;

road_bound_showFlag = true;  % 仅用于画图检查
doPlotScene = true;          % 是否先把场景画出来
scene = get_scene_from_json_withRoadBound(jsonPath, sampleIdx, road_bound_showFlag, doPlotScene);

lanes = scene.lane_centerlines;     % cell
road_bd = scene.road_boundaries;    % struct/cell/numeric（取决于你 JSON）
nv_past = scene.nv_past_xy;         % [Nh x 2]
nv_fut  = scene.nv_future_xy;       % [Nf x 2] (可为空)

% 组合 NV “全轨迹”（用于障碍车辆真实轨迹）
nv_all = nv_past;
if ~isempty(nv_fut)
    nv_all = [nv_all; nv_fut];
end

%% ========== 2) 虚构 AV 初始状态（相对 NV 偏移一点点） ==========
% NV 当前点（local frame）
p_nv0 = nv_past(end,:);

% NV 航向：用最后两个历史点估计；如果历史点不足，就用 0
if size(nv_past,1) >= 2
    dp = nv_past(end,:) - nv_past(end-1,:);
    theta0 = atan2(dp(2), dp(1));
else
    theta0 = 0;
end

% NV 速度：如果 JSON 里 past_traj 有 v，可从 m.past_traj(:,3) 取；否则用有限差分估计
v0 = estimate_speed_from_past(scene.sample, nv_past);

% AV 初始相对偏移（你可以改：横向 1.0~2.0m 常见）
offset_xy = [0.0, -1.5];  % [dx, dy] in local frame
p_av0 = p_nv0 + offset_xy;

% AV 初始状态 [x;y;v;theta]
S0 = [p_av0(1); p_av0(2); v0; theta0];

%% ========== 3) 参考轨迹：选最近的 lane centerline ==========
lane_ref = pick_nearest_lane(lanes, p_nv0);
% 将参考 lane 重采样成“按弧长均匀”的点序列，方便滚动取 N+1 点
ds = 0.5;  % 参考点间距（m），可改 0.2~1.0
wpt = resample_polyline_by_arclen(lane_ref, ds);

%% ========== 4) iLQR/滚动时域参数（简化版，不含 Σ） ==========
T_s = 0.1;
N   = 40;    % 4s horizon
T_sim = 8.0; % 你现在有 NV 8s 轨迹，先跑 8s 看看
K_exec = min(round(T_sim/T_s), size(nv_all,1)-1); % 确保 NV 轨迹够长

N_iter  = 8;
N_stage = 2;
lambda_damp = 1e-3;

% 车辆长度（双圆模型）
L_AV = 4.5;
L_NV = 4.5;

% 权重（你原来的）
w1 = 2.0; w2 = 0.1; w3 = 1.0; w4 = 3.0;

% 输入约束
a_min = -5; a_max = 3;
delta_lim = 0.5;
delta_min = -delta_lim; delta_max = delta_lim;

% 不膨胀安全距离：只用固定安全距离
s_safe = 2.4;
margin = 1.0;
s_safe_eff = s_safe + margin;

% 速度非负（保留）
v_min = 0.0;
eps_v = 1e-3;

% warm start
U_warm = zeros(2, N);

% 存储
S_exec = zeros(4, K_exec+1);
U_exec = zeros(2, K_exec);
S_exec(:,1) = S0;

%% ========== 5) 滚动时域求解 ==========
for k0 = 1:K_exec
    S_init = S_exec(:,k0);

    % --- 5.1) 取当前时刻 NV 的真实状态序列（horizon 内） ---
    % NV 直接用 nv_all(k,:) 的位置，速度/航向简单估计（够用来当动态障碍）
    nv_seg = build_nv_segment(nv_all, k0, N);

    % --- 5.2) 从 lane wpt 中截取 N+1 个参考点（跟随“道路中心线”） ---
    % 用“离当前 AV 最近”的参考点索引做起点（避免永远从头开始）
    wpt_seg = slice_wpt_for_horizon(wpt, [S_init(1), S_init(2)], N);

    x_wpt = wpt_seg(:,1);
    y_wpt = wpt_seg(:,2);

    % --- 5.3) 求解 HU-CILQR（这里是：双圆避障 + road boundary barrier） ---
    [S_pred, U_opt] = solve_HU_CILQR_noSigma_withRoadBound( ...
        S_init, U_warm, x_wpt, y_wpt, v0, ...
        nv_seg, L_AV, L_NV, s_safe_eff, ...
        road_bd, ...
        w1,w2,w3,w4, a_min,a_max, delta_min,delta_max, v_min, eps_v, ...
        T_s, N, N_stage, N_iter, lambda_damp);

    % --- 5.4) 执行第一个控制 ---
    U0 = U_opt(:,1);
    U_exec(:,k0) = U0;
    S_exec(:,k0+1) = f_dynamics(S_exec(:,k0), U0, T_s);

    U_warm = [U_opt(:,2:end), U_opt(:,end)];
end

%% ========== 6) 结果绘图（论文风格） ==========
figure('Name','AV plan vs NV GT + Lanes + RoadBound'); hold on; grid on; axis equal;

% lanes（灰）
for l = 1:numel(lanes)
    xy = lanes{l};
    if isempty(xy), continue; end
    plot(xy(:,1), xy(:,2), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
end

% road boundaries（深灰，可开关）
if road_bound_showFlag && ~isempty(road_bd)
    plot_road_boundaries_local(road_bd);
end

% NV（红虚线）
plot(nv_all(1:K_exec+1,1), nv_all(1:K_exec+1,2), 'r--', 'LineWidth', 2.0);

% AV（蓝实线）
plot(S_exec(1,:), S_exec(2,:), 'b-', 'LineWidth', 2.5);

xlabel('x_{rel} (m)'); ylabel('y_{rel} (m)');
title(sprintf('Sample #%d: AV (planned) vs NV (GT)', sampleIdx));

legend({'Lane centerlines','Road boundaries','NV (GT)','AV (planned)'}, 'Location','bestoutside');

%% ========================= 辅助函数 =========================

function v0 = estimate_speed_from_past(m, nv_past_xy)
% 优先从 past_traj 第3列 v 取；否则差分估计
v0 = 5.0;
try
    if isfield(m,'past_traj') && ~isempty(m.past_traj) && size(m.past_traj,2) >= 3
        v0 = m.past_traj(end,3);
    else
        if size(nv_past_xy,1) >= 2
            dp = nv_past_xy(end,:) - nv_past_xy(end-1,:);
            v0 = norm(dp) / 0.5; % nuScenes 2Hz => 0.5s
        end
    end
catch
end
end

function lane_ref = pick_nearest_lane(lanes, p0)
% 选取离 p0 最近的 lane polyline
bestD = inf; best = [];
for i = 1:numel(lanes)
    xy = lanes{i};
    if isempty(xy), continue; end
    d = min(vecnorm(xy - p0, 2, 2));
    if d < bestD
        bestD = d;
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

function nv_seg = build_nv_segment(nv_all, k0, N)
% NV 轨迹段：用位置 (x,y)，速度/航向简单差分补齐
idx0 = k0;
idx1 = min(k0+N, size(nv_all,1));
xy = nv_all(idx0:idx1,:);
if size(xy,1) < N+1
    xy = [xy; repmat(xy(end,:), N+1-size(xy,1), 1)];
end

v = zeros(N,1);
th = zeros(N,1);
for k = 1:N
    dp = xy(k+1,:) - xy(k,:);
    th(k) = atan2(dp(2), dp(1));
    v(k)  = norm(dp)/0.1; % 这里按 T_s=0.1 近似（够做动态障碍）
end
nv_seg = struct('xy',xy,'v',v,'th',th);
end

function S_next = f_dynamics(S, U, T_s)
S_next = S + [S(3)*cos(S(4));
              S(3)*sin(S(4));
              U(1);
              U(2)] * T_s;
end

function [S_traj, U_best] = solve_HU_CILQR_noSigma_withRoadBound( ...
    S0, U_init, x_wpt, y_wpt, v_d, ...
    nv_seg, L_AV, L_NV, s_safe_eff, ...
    road_bd, ...
    w1,w2,w3,w4, a_min,a_max, delta_min,delta_max, v_min, eps_v, ...
    T_s, N, N_stage, N_iter, lambda_damp)

% 这里为了“能跑 + 对齐结构”，给一个简化 iLQR：
% - running cost: w1*(pos-wpt)^2 + w2*(v-vd)^2 + w3*a^2 + w4*delta^2
% - obstacle barrier: exp(-(dmin - s_safe_eff))
% - road boundary barrier: outside polygon => 大罚；inside near boundary => exp(-dist)
%
% 说明：这段你之后可以直接替换成你 HU_CILQR_Simple 里的 solve 函数框架。
U = U_init;
S_traj = zeros(4,N+1);

% 预构造 road polygon（如果能构造出来）
roadPoly = build_road_polyshape(road_bd);

for stage = 1:N_stage
    for it = 1:N_iter
        % forward rollout
        S_traj(:,1) = S0;
        for k = 1:N
            S_traj(:,k+1) = f_dynamics(S_traj(:,k), U(:,k), T_s);
        end

        % backward pass（极简近似：只做一步比例下降，避免写满你整套 iLQR）
        % 你如果希望保留完整 iLQR，我建议直接把 HU_CILQR_Simple 的 solve 函数复制过来，
        % 然后把 Σ / s_safe_eff 的部分删掉，并把 road barrier 加进 running_cost 即可。
        alpha = 0.3;
        for k = 1:N
            % 用数值梯度做一个“够用”的下降方向（2维输入）
            u0 = U(:,k);
            g = numeric_grad_u(@(uu) running_cost_one(S_traj(:,k), uu, ...
                    [x_wpt(k);y_wpt(k)], v_d, ...
                    nv_seg.xy(k,:), L_AV, L_NV, s_safe_eff, ...
                    roadPoly, stage, v_min, eps_v, w1,w2,w3,w4), u0);

            U(:,k) = U(:,k) - alpha*g;

            % clamp
            U(1,k) = min(max(U(1,k), a_min), a_max);
            U(2,k) = min(max(U(2,k), delta_min), delta_max);
        end
    end
end

U_best = U;
end

function J = running_cost_one(S, U, wpt, v_d, nv_xy, L_AV, L_NV, s_safe_eff, roadPoly, stage, v_min, eps_v, w1,w2,w3,w4)
% base quadratic
ex = S(1)-wpt(1);
ey = S(2)-wpt(2);
ev = S(3)-v_d;
J = w1*(ex^2 + ey^2) + w2*(ev^2) + w3*(U(1)^2) + w4*(U(2)^2);

% v>=0 barrier
if stage==1, wv=10; else, wv=80; end
den = (S(3)-v_min+eps_v);
if den < 0.5
    J = J + wv/den;
end

% obstacle barrier (NV as moving obstacle): double-circle min dist
dmin = doublecircle_dmin_xytheta(S, nv_xy, L_AV, L_NV);
gap = dmin - s_safe_eff;
if stage==1, wobs=30; range=6; else, wobs=150; range=8; end
if gap < range
    J = J + wobs*exp(-gap);
end

% road boundary barrier
if ~isempty(roadPoly)
    [in,on] = isinterior(roadPoly, S(1), S(2));
    if ~(in || on)
        % outside: heavy penalty
        if stage==1, wb=200; else, wb=2000; end
        J = J + wb*(1 + min_distance_to_polyshape_boundary(roadPoly, [S(1),S(2)]));
    else
        % inside but near boundary: soft repulsion
        d2bd = min_distance_to_polyshape_boundary(roadPoly, [S(1),S(2)]);
        if stage==1, wb=5; else, wb=30; end
        if d2bd < 2.0
            J = J + wb*exp(-(d2bd)); % 越靠近边界惩罚越大
        end
    end
end
end

function g = numeric_grad_u(fun, u)
% 数值梯度（2维）
eps = 1e-3;
g = zeros(2,1);
for i=1:2
    du = zeros(2,1); du(i)=eps;
    g(i) = (fun(u+du)-fun(u-du))/(2*eps);
end
end

function dmin = doublecircle_dmin_xytheta(S_av, xy_nv, L_AV, L_NV)
% 用 AV 的 (x,y,theta) 与 NV 的 (x,y) 计算双圆最小距离（NV theta 不用也行）
x = S_av(1); y=S_av(2); th=S_av(4);
Fe = [x;y] + (L_AV/2)*[cos(th); sin(th)];
Be = [x;y] - (L_AV/2)*[cos(th); sin(th)];

% NV：用“等效朝向=0”也可以（你有 nv_seg.th 的话也能加进去）
Fo = xy_nv(:) + (L_NV/2)*[1;0];
Bo = xy_nv(:) - (L_NV/2)*[1;0];

dmin = min([norm(Fe-Fo), norm(Fe-Bo), norm(Be-Fo), norm(Be-Bo)]);
end

function roadPoly = build_road_polyshape(road_bd)
% 尝试从 road_boundaries 构造一个 polyshape（假设是左右边界两条折线）
roadPoly = [];
try
    if iscell(road_bd) && numel(road_bd) >= 2
        L = road_bd{1}; R = road_bd{2};
        if size(L,2)==2 && size(R,2)==2
            % 拼成闭合多边形：L 正向 + R 反向
            P = [L; flipud(R)];
            roadPoly = polyshape(P(:,1), P(:,2));
        end
    elseif isstruct(road_bd) && numel(road_bd) >= 2
        % 你也可能存成 struct(points=[])
        L = road_bd(1).points; R = road_bd(2).points;
        P = [L; flipud(R)];
        roadPoly = polyshape(P(:,1), P(:,2));
    end
catch
    roadPoly = [];
end
end

function d = min_distance_to_polyshape_boundary(poly, p)
% 计算点到 polyshape 边界的最小距离（近似：采样边界点）
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
