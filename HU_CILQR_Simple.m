%% HU-CILQR-lite（双圆车-车）+ 滚动时域规划（0~18s）
% 目标：提供“结构对齐”的可运行骨架，便于后续接入：
%   (1) 车道线约束  (2) 真实 NV 轨迹  (3) Cantelli 机会约束/均值方差聚合
%
% --- 规划设置（对应论文/Markdown）---
% 规划步长：T_s = 0.1s
% 规划时域：N = 40 steps  =>  T_h = N*T_s = 4s
% 滚动规划：0~18s，总执行步 K_exec = 18/0.1 = 180
%
% --- 模型（对应 Markdown）---
% 状态 S_t = [x_t; y_t; v_t; theta_t]
% 输入 U_t = [a_t; delta_t]  (a:加速度, delta:角速度/yaw rate)
%
% --- 代价函数（对应你 Markdown 的二次型）---
% J = 终端 + 累积：
%   w1*||p_t - p_t^wpt||^2 + w2*(v_t - v^d)^2 + w3*a_t^2 + w4*delta_t^2
%
% --- 安全约束（当前为“工程近似版”）---
% 双圆模型 4 对圆心最小距离 d_min
% 安全距离膨胀：s_safe_eff(t) = s_safe + beta*sqrt(lambda_max(Sigma_t)) + margin
% clearance(t) = d_min(t) - s_safe_eff(t)
% （注意：这还不是 Cantelli 机会约束的 μ/σ² 形式，仅为近似替代项）

clc; clear; close all;

%% ========================= 1) 参数设置（Markdown 符号对齐）=========================
T_s   = 0.1;          % 采样周期/时间步长  T_s
N     = 40;           % 规划时域步数        N
T_sim = 18.0;         % 总仿真时长          0~18s
K_exec = round(T_sim / T_s);     % 实际执行步数（滚动规划次数）

N_iter  = 8;          % 每阶段 iLQR 迭代次数
N_stage = 2;          % 1:软约束 2:硬约束（barrier 权重更大）
lambda_damp = 1e-3;   % iLQR 阻尼（Quu 正则）

% 初始状态 S_{t_sta}
S0 = [0; 0; 5; 0];    % [x0;y0;v0;theta0]

% 车辆长度（双圆模型：前/后圆心 = +/- L/2）
L_AV = 4.5;           % AV 车长
L_NV = 4.5;           % NV 车长（示例）

% 目标速度 v^d
v_d = 9.0;

% 权重（对应 Markdown 的两种 case）
w1 = 2.0;
w2 = 0.1;   % NV 换道影响 AV：0.1；若 AV 主动换道可取 0.8
w3 = 1.0;
w4 = 3.0;

% 输入约束（Markdown 里你设定的典型范围）
a_min = -5;  a_max = 3;          % a_t ∈ [a_min, a_max]
delta_lim = 0.5;                 % |delta_t| ≤ delta_lim
delta_min = -delta_lim; delta_max = delta_lim;

% 安全距离参数（Markdown）
s_safe = 2.4;        % s_safe
beta   = 2.0;        % beta
margin = 1.0;        % margin

% 速度非负（工程补丁：避免倒车）
v_min = 0.0;
eps_v = 1e-3;

%% ========================= 2) 参考轨迹 p_t^{wpt}（直线示例）=========================
% 预生成：整个仿真 + 额外 N 步（避免滚动切片越界）
x_wpt_all = linspace(S0(1), S0(1) + v_d*(T_sim + N*T_s), K_exec + N + 1);
y_wpt_all = zeros(size(x_wpt_all));

%% ========================= 3) NV 状态均值 \hat{S}_t^i + 协方差 Σ_t（示例）===========
% 注意：这里仍是“人造 NV”，后续替换成 nuScenes/HighD/INTERACTION 的真实轨迹即可
S_hat_NV_all = zeros(4, K_exec + N);  % \hat{S}_t = [\hat{x};\hat{y};\hat{v};\hat{theta}]
v_NV = 2.0;
y_NV = 0.5;
x_NV0 = 15.0;

for k = 1:(K_exec + N)
    S_hat_NV_all(:,k) = [x_NV0 + v_NV*(k-1)*T_s;
                         y_NV;
                         v_NV;
                         0.0];
end

% Σ_t：仅对位置 (x,y) 设 2x2（示例随时间增大）
Sigma_all = zeros(2,2, K_exec + N);
for k = 1:(K_exec + N)
    sx = 0.20 + 0.02*(k-1);
    sy = 0.10 + 0.01*(k-1);
    Sigma_all(:,:,k) = diag([sx^2, sy^2]);
end

%% ========================= 4) 滚动规划存储（执行轨迹）=========================
S_exec = zeros(4, K_exec+1);
U_exec = zeros(2, K_exec);
S_exec(:,1) = S0;

% warm-start：每次滚动规划的控制序列 U_{t_sta:t_plan-1}
U_warm = zeros(2, N);

% 日志：每次滚动规划得到的最优代价 + horizon 内最小安全裕度
J_roll = zeros(1, K_exec);
minClear_roll = zeros(1, K_exec);

%% ========================= 5) 滚动时域主循环（Receding Horizon）====================
for k0 = 1:K_exec
    % 当前时刻的初始状态：S_{t_sta}
    S_init = S_exec(:,k0);

    % 取出 horizon 内的参考点：p_t^{wpt}
    idx_wpt = k0:(k0+N);  % N+1 点（含终端）
    x_wpt = x_wpt_all(idx_wpt);
    y_wpt = y_wpt_all(idx_wpt);

    % 取出 horizon 内的 NV 均值轨迹与协方差：\hat{S}_t, Σ_t
    idx_h = k0:(k0+N-1);  % N 步（过程代价/约束）
    S_hat_NV = S_hat_NV_all(:, idx_h);
    Sigma_h  = Sigma_all(:,:, idx_h);

    % HU-CILQR-lite 求解（双圆 + 安全距离膨胀 + barrier）
    [S_pred, U_opt, J_best, minClear] = solve_HU_CILQR_lite( ...
        S_init, U_warm, ...
        x_wpt, y_wpt, v_d, ...
        S_hat_NV, Sigma_h, ...
        L_AV, L_NV, ...
        s_safe, beta, margin, ...
        w1,w2,w3,w4, ...
        a_min,a_max, delta_min,delta_max, ...
        v_min, eps_v, ...
        T_s, N, N_stage, N_iter, lambda_damp);

    % 只执行第一个控制（滚动规划核心）
    U0 = U_opt(:,1);
    U_exec(:,k0) = U0;

    % 状态推进一步：S_{t+1} = f(S_t, U_t)
    S_exec(:,k0+1) = f_dynamics(S_exec(:,k0), U0, T_s);

    % warm-start 移位：丢掉第 1 个，末尾复制最后一个
    U_warm = [U_opt(:,2:end), U_opt(:,end)];

    % log
    J_roll(k0) = J_best;
    minClear_roll(k0) = minClear;

    if mod(k0,20)==0
        fprintf('[%3d/%3d] J=%.2f, minClear(horizon)=%.3f\n', k0, K_exec, J_best, minClear);
    end
end

%% ========================= 6) 评估与绘图（以符号解释变量）=========================
t   = (0:K_exec)*T_s;
t_u = (0:K_exec-1)*T_s;

x = S_exec(1,:);  y = S_exec(2,:);  v = S_exec(3,:);  theta = S_exec(4,:);
a = U_exec(1,:);  delta = U_exec(2,:);

% 计算每个执行时刻的最小距离 d_min 与 s_safe_eff
d_min_ts = zeros(1, K_exec);
s_safe_eff_ts = zeros(1, K_exec);

for k = 1:K_exec
    S_hat = S_hat_NV_all(:,k);
    Sig   = Sigma_all(:,:,k);

    sigma_max = sqrt(max(eig(Sig)));
    s_safe_eff = s_safe + beta*sigma_max + margin;

    d_min_ts(k) = doublecircle_dmin(S_exec(:,k), S_hat, L_AV, L_NV);
    s_safe_eff_ts(k) = s_safe_eff;
end

clearance = d_min_ts - s_safe_eff_ts;
[min_clear_all, idx_mc] = min(clearance);

% XY
figure('Name','XY: AV vs NV(mean)'); grid on; hold on; axis equal;
plot(x, y, 'b-', 'LineWidth', 2);
plot(S_hat_NV_all(1,1:K_exec), S_hat_NV_all(2,1:K_exec), 'r--', 'LineWidth', 2);
title('XY 轨迹：AV 与 NV（均值）');
xlabel('x (m)'); ylabel('y (m)');
legend('AV','NV mean','Location','best');

% 画几个膨胀安全圈（仅示意）
show_idx = unique(round(linspace(1, K_exec, 5)));
for kk = show_idx
    Sig = Sigma_all(:,:,kk);
    sigma_max = sqrt(max(eig(Sig)));
    s_safe_eff = s_safe + beta*sigma_max + margin;
    viscircles(S_hat_NV_all(1:2,kk)', s_safe_eff, 'Color',[1 0 0], 'LineStyle',':');
end

% 评估面板
figure('Name','Evaluation Panel (5x1)');
subplot(5,1,1); grid on; plot(t, v, 'LineWidth',1.4);
ylabel('v_t (m/s)'); xlim([t(1), t(end)]);

subplot(5,1,2); grid on; hold on;
plot(t_u, d_min_ts, 'LineWidth',1.4);
plot(t_u, s_safe_eff_ts, '--', 'LineWidth',1.4);
ylabel('d_{min} / s_{safe,eff}'); xlim([t_u(1), t_u(end)]);
legend('d_{min}(t)','s_{safe,eff}(t)','Location','best');

subplot(5,1,3); grid on; plot(t_u, a, 'LineWidth',1.4);
ylabel('a_t (m/s^2)'); xlim([t_u(1), t_u(end)]);

subplot(5,1,4); grid on; plot(t_u, delta, 'LineWidth',1.4);
ylabel('\delta_t (rad/s)'); xlim([t_u(1), t_u(end)]);

subplot(5,1,5); grid on; hold on;
plot(t_u, clearance, 'LineWidth',1.4);
plot(t_u(idx_mc), min_clear_all, 'o', 'MarkerSize',7, 'LineWidth',2);
yline(0,'--','LineWidth',1.1);
ylabel('clearance'); xlabel('t (s)');
xlim([t_u(1), t_u(end)]);
legend('d_{min}-s_{safe,eff}','min','0','Location','best');

figure('Name','Rolling Best Cost'); grid on;
plot(t_u, J_roll, 'LineWidth',1.4);
xlabel('t (s)'); ylabel('J'); title('每次滚动规划的最优总代价');

fprintf('\n=== Summary ===\n');
fprintf('Global min clearance = %.3f m at t=%.2f s\n', min_clear_all, t_u(idx_mc));
fprintf('Final: (x,y)=(%.2f,%.2f), v=%.2f\n', x(end), y(end), v(end));

%% ========================= 子函数（HU-CILQR-lite）==========================

function S_next = f_dynamics(S, U, T_s)
% 离散动力学：S_{t+1} = f(S_t, U_t)
% S=[x;y;v;theta], U=[a;delta]
    S_next = S + [S(3)*cos(S(4));
                  S(3)*sin(S(4));
                  U(1);
                  U(2)] * T_s;
end

function [A,B] = linearize_f(S, ~, T_s)
% 一阶泰勒线性化：δS_{t+1} ≈ A_t δS_t + B_t δU_t
    A = [1, 0, T_s*cos(S(4)), -T_s*S(3)*sin(S(4));
         0, 1, T_s*sin(S(4)),  T_s*S(3)*cos(S(4));
         0, 0, 1, 0;
         0, 0, 0, 1];
    B = [0, 0;
         0, 0;
         T_s, 0;
         0, T_s];
end

function [S_traj, U_best, Jbest, minClear] = solve_HU_CILQR_lite( ...
    S0, U_init, x_wpt, y_wpt, v_d, ...
    S_hat_NV, Sigma_h, ...
    L_AV, L_NV, ...
    s_safe, beta, margin, ...
    w1,w2,w3,w4, ...
    a_min,a_max, delta_min,delta_max, ...
    v_min, eps_v, ...
    T_s, N, N_stage, N_iter, lambda_damp)

% HU-CILQR-lite：
% - 两阶段 barrier（软->硬）
% - 双圆车-车最小距离 d_min
% - s_safe_eff(t) 随 Σ_t 膨胀（工程近似）
% 输出：
%   S_traj: horizon 轨迹
%   U_best: 最优控制序列
%   Jbest: 最优总代价
%   minClear: horizon 内最小 clearance

    U = U_init;
    S_traj = zeros(4, N+1);

    Jbest = inf;
    U_best = U;
    minClear = inf;

    for stage = 1:N_stage
        for it = 1:N_iter
            % 前向 rollout
            S_traj(:,1) = S0;
            for k = 1:N
                S_traj(:,k+1) = f_dynamics(S_traj(:,k), U(:,k), T_s);
            end

            % 计算总代价 + 最小裕度
            [J, minClrNow] = total_cost_and_clearance( ...
                S_traj, U, x_wpt, y_wpt, v_d, S_hat_NV, Sigma_h, ...
                L_AV, L_NV, s_safe, beta, margin, ...
                w1,w2,w3,w4, stage, v_min, eps_v);

            if J < Jbest
                Jbest = J;
                U_best = U;
                minClear = minClrNow;
            end

            % 终端导数
            [Vx, Vxx] = terminal_derivatives(S_traj(:,N+1), x_wpt(end), y_wpt(end), v_d, w1, w2);

            % iLQR backward
            K = zeros(2,4,N);
            kff = zeros(2,N);

            for k = N:-1:1
                Sk = S_traj(:,k);
                Uk = U(:,k);

                [A,B] = linearize_f(Sk, Uk, T_s);

                [lx,lxx,lu,luu,lux] = running_derivatives( ...
                    Sk, Uk, x_wpt(k), y_wpt(k), v_d, ...
                    S_hat_NV(:,k), Sigma_h(:,:,k), ...
                    L_AV, L_NV, s_safe, beta, margin, ...
                    w1,w2,w3,w4, stage, v_min, eps_v);

                Qx  = lx  + A' * Vx;
                Qu  = lu  + B' * Vx;
                Qxx = lxx + A' * Vxx * A;
                Quu = luu + B' * Vxx * B;
                Qux = lux + B' * Vxx * A;

                Quu_reg = Quu + lambda_damp*eye(2);

                kff(:,k) = -(Quu_reg \ Qu);
                K(:,:,k) = -(Quu_reg \ Qux);

                Vx  = Qx  + K(:,:,k)'*Quu*kff(:,k) + K(:,:,k)'*Qu + Qux'*kff(:,k);
                Vxx = Qxx + K(:,:,k)'*Quu*K(:,:,k) + K(:,:,k)'*Qux + Qux'*K(:,:,k);
                Vxx = 0.5*(Vxx + Vxx');
            end

            % 前馈更新（固定 alpha）
            alpha = 0.5;
            dS = zeros(4,1);

            S_new = zeros(4,N+1);
            S_new(:,1) = S0;
            U_new = U;

            for k = 1:N
                dU = kff(:,k) + K(:,:,k)*dS;
                Uk = U(:,k) + alpha*dU;

                % 输入裁剪：a_t、delta_t
                Uk(1) = min(max(Uk(1), a_min), a_max);
                Uk(2) = min(max(Uk(2), delta_min), delta_max);

                S_next = f_dynamics(S_new(:,k), Uk, T_s);

                dS = S_next - S_traj(:,k+1);

                U_new(:,k) = Uk;
                S_new(:,k+1) = S_next;
            end

            U = U_new;
            S_traj = S_new;
        end
    end

    % 返回 best
    U = U_best;
    S_traj(:,1) = S0;
    for k = 1:N
        S_traj(:,k+1) = f_dynamics(S_traj(:,k), U(:,k), T_s);
    end
end

function [Vx, Vxx] = terminal_derivatives(S_T, xw_T, yw_T, v_d, w1, w2)
% 终端代价：w1||p_T - p_T^{wpt}||^2 + w2(v_T - v^d)^2
    ex = S_T(1)-xw_T;
    ey = S_T(2)-yw_T;
    ev = S_T(3)-v_d;

    Vx  = [2*w1*ex; 2*w1*ey; 2*w2*ev; 0];
    Vxx = diag([2*w1, 2*w1, 2*w2, 0]);
end

function [lx,lxx,lu,luu,lux] = running_derivatives( ...
    S, U, xw, yw, v_d, S_hat, Sigma, ...
    L_AV, L_NV, s_safe, beta, margin, ...
    w1,w2,w3,w4, stage, v_min, eps_v)

% 过程代价（Markdown 二次型） + barrier（工程近似）
% 1) 二次型：w1||p - p^{wpt}||^2 + w2(v-v^d)^2 + w3 a^2 + w4 delta^2
% 2) v>=0 barrier：防止倒车（不是论文重点，但工程上必要）
% 3) 双圆安全距离 barrier：用 s_safe_eff 膨胀（近似替代机会约束）

    ex = S(1)-xw;  ey = S(2)-yw;  ev = S(3)-v_d;

    lx  = [2*w1*ex; 2*w1*ey; 2*w2*ev; 0];
    lxx = diag([2*w1, 2*w1, 2*w2, 0]);

    lu  = [2*w3*U(1); 2*w4*U(2)];
    luu = diag([2*w3, 2*w4]);
    lux = zeros(2,4);

    % ---------- v>=0 barrier ----------
    if stage==1, w_vbar = 10; else, w_vbar = 80; end
    denom = (S(3) - v_min + eps_v);
    if denom < 0.5
        lx(3) = lx(3) - w_vbar/(denom^2);
        lxx(3,3) = lxx(3,3) + 2*w_vbar/(denom^3);
    end

    % ---------- 双圆安全距离 barrier（膨胀） ----------
    sigma_max = sqrt(max(eig(Sigma)));
    s_safe_eff = s_safe + beta*sigma_max + margin;

    if stage==1
        w_obs = 20; sense = 5.0;
    else
        w_obs = 120; sense = 7.0;
    end

    [dmin, grad_xyz] = dmin_grad_xytheta(S, S_hat, L_AV, L_NV); % 对 [x;y;theta] 的梯度
    gap = dmin - s_safe_eff;

    if gap < sense
        b = w_obs * exp(-(gap));
        lx(1) = lx(1) - b * grad_xyz(1);
        lx(2) = lx(2) - b * grad_xyz(2);
        lx(4) = lx(4) - b * grad_xyz(3);

        g = [grad_xyz(1); grad_xyz(2); 0; grad_xyz(3)];
        lxx = lxx + b*(g*g');  % Hessian 近似
        lxx = 0.5*(lxx + lxx');
    end
end

function [J, minClear] = total_cost_and_clearance( ...
    S_traj, U, x_wpt, y_wpt, v_d, S_hat_NV, Sigma_h, ...
    L_AV, L_NV, s_safe, beta, margin, w1,w2,w3,w4, stage, v_min, eps_v)

% 计算总代价 J，并记录 horizon 内最小安全裕度
    N = size(U,2);
    J = 0;
    minClear = inf;

    for k = 1:N
        S = S_traj(:,k);
        Uk = U(:,k);

        ex = S(1)-x_wpt(k); ey = S(2)-y_wpt(k); ev = S(3)-v_d;
        J = J + w1*(ex^2+ey^2) + w2*(ev^2) + w3*(Uk(1)^2) + w4*(Uk(2)^2);

        if stage==1, w_vbar=10; w_obs=20; sense=5.0;
        else,        w_vbar=80; w_obs=120; sense=7.0;
        end

        denom = (S(3) - v_min + eps_v);
        if denom < 0.5
            J = J + w_vbar/denom;
        end

        Sigma = Sigma_h(:,:,k);
        sigma_max = sqrt(max(eig(Sigma)));
        s_safe_eff = s_safe + beta*sigma_max + margin;

        dmin = doublecircle_dmin(S, S_hat_NV(:,k), L_AV, L_NV);
        gap = dmin - s_safe_eff;
        minClear = min(minClear, gap);

        if gap < sense
            J = J + w_obs*exp(-(gap));
        end
    end

    % terminal
    S_T = S_traj(:,end);
    exT = S_T(1)-x_wpt(end); eyT = S_T(2)-y_wpt(end); evT = S_T(3)-v_d;
    J = J + w1*(exT^2+eyT^2) + w2*(evT^2);
end

function dmin = doublecircle_dmin(S_AV, S_NV, L_AV, L_NV)
% 双圆模型的 4 对圆心距离最小值 d_min
    [F_AV, B_AV] = front_back_centers(S_AV, L_AV);
    [F_NV, B_NV] = front_back_centers(S_NV, L_NV);

    dFF = norm(F_AV - F_NV);
    dFB = norm(F_AV - B_NV);
    dBF = norm(B_AV - F_NV);
    dBB = norm(B_AV - B_NV);

    dmin = min([dFF, dFB, dBF, dBB]);
end

function [dmin, grad_xyz] = dmin_grad_xytheta(S_AV, S_NV, L_AV, L_NV)
% d_min 对 AV [x;y;theta] 的梯度（分段：对激活的最小距离那对求导）
    [F_AV, B_AV, dF_dth, dB_dth] = front_back_centers_with_jac(S_AV, L_AV);
    [F_NV, B_NV] = front_back_centers(S_NV, L_NV);

    P = {F_AV,F_NV,dF_dth;  F_AV,B_NV,dF_dth;  B_AV,F_NV,dB_dth;  B_AV,B_NV,dB_dth};
    d = zeros(1,4);
    for i=1:4
        d(i) = norm(P{i,1} - P{i,2});
    end
    [dmin, idx] = min(d);

    p = P{idx,1}; q = P{idx,2}; dp_dth = P{idx,3};
    r = p - q;
    dist = max(norm(r), 1e-6);
    grad_p = r/dist;                 % ∂dist/∂p

    dd_dx = grad_p(1);
    dd_dy = grad_p(2);
    dd_dth = grad_p' * dp_dth;

    grad_xyz = [dd_dx; dd_dy; dd_dth];
end

function [F,B] = front_back_centers(S, L)
% 前/后圆心：p ± (L/2)[cosθ; sinθ]
    x = S(1); y = S(2); th = S(4);
    F = [x;y] + (L/2)*[cos(th); sin(th)];
    B = [x;y] - (L/2)*[cos(th); sin(th)];
end

function [F,B,dF_dth,dB_dth] = front_back_centers_with_jac(S, L)
% 前/后圆心 + 对 θ 的导数
    x = S(1); y = S(2); th = S(4);
    c = cos(th); s = sin(th);

    F = [x;y] + (L/2)*[c; s];
    B = [x;y] - (L/2)*[c; s];

    dF_dth = (L/2)*[-s; c];
    dB_dth = -(L/2)*[-s; c];
end
