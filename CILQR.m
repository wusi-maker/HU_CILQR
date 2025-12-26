%% HU-CILQR-lite in MATLAB (Enhanced: 3D traj + subplot panel + risk & convergence)
% - Dynamic obstacle mean
% - Time-varying obstacle covariance Sigma(k)
% - Distance inflation: d_safe_eff(k) = d_safe + beta * sqrt(lambda_max(Sigma(k))) + margin
% - Two-stage: soft -> hard
% - Damped Quu
% - Velocity non-negative barrier
% - Outputs: traj3d = [x, y, t]
% - Plots:
%   Fig1: XY path
%   Fig2: 3D (x,y,t)
%   Fig3: subplot(5,1): v(t), dist & d_safe_eff, a(t), omega(t), clearance(t) + cost trace

clc; clear; close all;

%% 1) Parameters
N        = 50;
dt       = 0.2;
max_iter = 10;
stages   = 2;

x0    = [0; 0; 5; 0];
u_nom = zeros(2, N);

v_target = 5.0;
y_target = 0.0;
x_goal   = x0(1) + v_target * N * dt;

d_safe = 2.4;
beta   = 2.0;
margin = 1.0;

lambda_damp = 1e-3;

a_min = -3.0; a_max = 3.0;
w_min = -0.8; w_max = 0.8;

% --- weights used for cost logging (must match running_cost) ---
Q_run = diag([0.01, 0.5, 0.5, 0.1]);
R_run = diag([0.1, 1.0]);
Qf    = diag([20, 20, 1, 1]);
w_progress = 0.1;

%% 2) Dynamic obstacle mean
obs_traj = zeros(2, N);
v_obs    = 2.0;
y_obs    = 0.5;
for k = 1:N
    obs_traj(:,k) = [15 + v_obs * (k-1) * dt; y_obs];
end

%% 3) Obstacle covariance sequence
Sigma_obs = zeros(2,2,N);
for k = 1:N
    sx = 0.2 + 0.05*(k-1);
    sy = 0.1 + 0.03*(k-1);
    Sigma_obs(:,:,k) = diag([sx^2, sy^2]);
end

%% 4) Storage
x_traj = zeros(4, N+1);
x_traj(:,1) = x0;

% For convergence logging
J_hist = [];  % append each (stage,iter) total cost
J_tag  = [];  % same length, stores stage number (1 or 2)

%% 5) Optimization (2-stage)
for stage = 1:stages
    for iter = 1:max_iter

        % --- Forward pass ---
        x_traj(:,1) = x0;
        for k = 1:N
            x_traj(:, k+1) = vehicle_dynamics(x_traj(:,k), u_nom(:,k), dt);
        end

        % --- Log total cost for this rollout (for convergence plot) ---
        J = compute_total_cost(x_traj, u_nom, obs_traj, Sigma_obs, ...
                               x_goal, y_target, v_target, ...
                               d_safe, beta, margin, stage, ...
                               Q_run, R_run, Qf, w_progress);
        J_hist(end+1) = J; %#ok<SAGROW>
        J_tag(end+1)  = stage;

        % --- Backward pass ---
        [lxT, lxxT] = terminal_cost(x_traj(:,N+1), x_goal, y_target, v_target, Qf);
        V_x  = lxT;
        V_xx = lxxT;

        K_gains = zeros(2, 4, N);
        k_ff    = zeros(2, N);

        for k = N:-1:1
            [Ak, Bk] = get_linearized_dynamics(x_traj(:,k), u_nom(:,k), dt);

            obs_k = obs_traj(:,k);
            Sig_k = Sigma_obs(:,:,k);

            [lx, lxx, lu, luu, lux] = running_cost( ...
                x_traj(:,k), u_nom(:,k), ...
                obs_k, Sig_k, ...
                x_goal, y_target, v_target, ...
                d_safe, beta, margin, stage, ...
                Q_run, R_run, w_progress);

            Qx  = lx  + Ak' * V_x;
            Qu  = lu  + Bk' * V_x;
            Qxx = lxx + Ak' * V_xx * Ak;
            Quu = luu + Bk' * V_xx * Bk;
            Qux = lux + Bk' * V_xx * Ak;

            Quu_reg = Quu + lambda_damp * eye(2);

            k_ff(:,k)      = - (Quu_reg \ Qu);
            K_gains(:,:,k) = - (Quu_reg \ Qux);

            V_x  = Qx  + K_gains(:,:,k)' * Quu * k_ff(:,k) + K_gains(:,:,k)' * Qu + Qux' * k_ff(:,k);
            V_xx = Qxx + K_gains(:,:,k)' * Quu * K_gains(:,:,k) + K_gains(:,:,k)' * Qux + Qux' * K_gains(:,:,k);
            V_xx = 0.5 * (V_xx + V_xx');
        end

        % --- Control update ---
        alpha = 0.5;
        dx = zeros(4,1);
        for k = 1:N
            du = k_ff(:,k) + K_gains(:,:,k) * dx;
            u_new = u_nom(:,k) + alpha * du;

            u_new(1) = min(max(u_new(1), a_min), a_max);
            u_new(2) = min(max(u_new(2), w_min), w_max);

            x_new = vehicle_dynamics(x_traj(:,k), u_new, dt);
            dx = x_new - x_traj(:,k+1);

            u_nom(:,k)    = u_new;
            x_traj(:,k+1) = x_new;
        end
    end
end

%% 6) Post-process time series
t   = (0:N) * dt;
t_u = (0:N-1) * dt;

px = x_traj(1,:);
py = x_traj(2,:);
v  = x_traj(3,:);
theta = x_traj(4,:);

a = u_nom(1,:);
omega = u_nom(2,:);

dist_to_obs = zeros(1, N+1);
d_safe_eff_series = zeros(1, N+1);

for k = 1:N
    obs_k = obs_traj(:,k);
    dist_to_obs(k) = norm([px(k); py(k)] - obs_k);

    sigma_max  = sqrt(max(eig(Sigma_obs(:,:,k))));
    d_safe_eff_series(k) = d_safe + beta * sigma_max + margin;
end
dist_to_obs(N+1) = norm([px(N+1); py(N+1)] - obs_traj(:,N));
d_safe_eff_series(N+1) = d_safe_eff_series(N);

clearance = dist_to_obs - d_safe_eff_series;
[min_clear, idx_min_clear] = min(clearance);

traj3d = [px(:), py(:), t(:)];
assignin('base','traj3d',traj3d);

%% 7) Figure 1: XY path
figure('Name','XY Trajectory'); hold on; grid on; axis equal;
plot(px, py, 'b-o', 'LineWidth', 2);
plot(obs_traj(1,:), obs_traj(2,:), 'r--', 'LineWidth', 2);

show_idx = unique(round(linspace(1, N, 4)));
for kk = show_idx
    viscircles(obs_traj(:,kk)', d_safe_eff_series(kk), 'Color', [1 0 0], 'LineStyle', ':');
end
title('HU-CILQR-lite: XY path');
xlabel('X (m)'); ylabel('Y (m)');
legend('Ego trajectory','Obstacle mean traj','Inflated safety (samples)','Location','best');

%% 8) Figure 2: 3D trajectory (x,y,t)
figure('Name','3D Trajectory (x,y,t)'); grid on; hold on;
plot3(px, py, t, 'b-', 'LineWidth', 2);
plot3(obs_traj(1,:), obs_traj(2,:), (0:N-1)*dt, 'r--', 'LineWidth', 2);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Time (s)');
title('3D Trajectory: (x, y, time)');
legend('Ego (x,y,t)','Obstacle mean (x,y,t)','Location','best');
view(45, 25);

%% 9) Figure 3: One-screen evaluation panel (subplot 5x1)
figure('Name','Evaluation Panel (5x1)');

% (1) speed
subplot(5,1,1); grid on; hold on;
plot(t, v, 'LineWidth', 1.8);
ylabel('v (m/s)');
title('Evaluation Panel: v(t), dist(t), a(t), \omega(t), clearance(t) + cost');
xlim([t(1), t(end)]);

% (2) distance & safety threshold
subplot(5,1,2); grid on; hold on;
plot(t, dist_to_obs, 'LineWidth', 1.8);
plot(t, d_safe_eff_series, '--', 'LineWidth', 1.8);
ylabel('dist (m)');
legend('dist(ego,obs)','d\_safe\_eff','Location','best');
xlim([t(1), t(end)]);

% (3) acceleration
subplot(5,1,3); grid on; hold on;
plot(t_u, a, 'LineWidth', 1.8);
ylabel('a (m/s^2)');
xlim([t_u(1), t_u(end)]);

% (4) yaw rate
subplot(5,1,4); grid on; hold on;
plot(t_u, omega, 'LineWidth', 1.8);
ylabel('\omega (rad/s)');
xlim([t_u(1), t_u(end)]);

% (5) clearance + convergence
subplot(5,1,5); grid on; hold on;
plot(t, clearance, 'LineWidth', 1.8);
plot(t(idx_min_clear), min_clear, 'o', 'MarkerSize', 7, 'LineWidth', 2);
yline(0,'--','LineWidth',1.2);
ylabel('clearance (m)');
xlabel('Time (s)');
xlim([t(1), t(end)]);
legend('dist - d\_safe\_eff','min clearance','0 line','Location','best');

% Add convergence inset (simple overlay on a second axis)
ax1 = gca;
ax2 = axes('Position', ax1.Position, 'Color', 'none');
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.XLim = [1, numel(J_hist)];
ax2.YLimMode = 'auto';
hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, 1:numel(J_hist), J_hist, 'k-', 'LineWidth', 1.2);
ylabel(ax2, 'J (total cost)');
xlabel(ax2, 'iter index (stage+iter)');
ax2.Box = 'off';

%% Print summary
fprintf('Min clearance = %.3f m at t = %.2f s (idx=%d)\n', min_clear, t(idx_min_clear), idx_min_clear);
fprintf('traj3d exported to workspace: size = [%d x %d]\n', size(traj3d,1), size(traj3d,2));

%% ===================== Helper Functions =====================

function x_next = vehicle_dynamics(x, u, dt)
    x_next = x + [x(3)*cos(x(4));
                  x(3)*sin(x(4));
                  u(1);
                  u(2)] * dt;
end

function [A, B] = get_linearized_dynamics(x, ~, dt)
    A = [1, 0, dt*cos(x(4)), -dt*x(3)*sin(x(4));
         0, 1, dt*sin(x(4)),  dt*x(3)*cos(x(4));
         0, 0, 1, 0;
         0, 0, 0, 1];
    B = [0, 0;
         0, 0;
         dt, 0;
         0, dt];
end

function [lx, lxx, lu, luu, lux] = running_cost( ...
    x, u, obs, Sigma, ...
    x_goal, y_target, v_target, ...
    d_safe, beta, margin, stage, ...
    Q, R, w_progress)

    dx = [x(1)-x_goal; x(2)-y_target; x(3)-v_target; wrapToPiLocal(x(4))];

    lx  = Q * dx;
    lxx = Q;
    lu  = R * u;
    luu = R;
    lux = zeros(2,4);

    % progress shaping
    lx(1) = lx(1) - w_progress;

    % negative speed penalty (backup)
    w_vneg = 50.0;
    if x(3) < 0
        lx(3)    = lx(3) + 2*w_vneg*x(3);
        lxx(3,3) = lxx(3,3) + 2*w_vneg;
    end

    % inflated safety distance
    sigma_max  = sqrt(max(eig(Sigma)));
    d_safe_eff = d_safe + beta*sigma_max + margin;

    r    = x(1:2) - obs;
    dist = max(norm(r), 1e-6);

    % stage params
    v = x(3);
    v_min = 0.05;
    epsv  = 1e-3;

    if stage == 1
        weight_obs  = 30;
        sense_range = 4.0;
        w_vbar      = 30;
    else
        weight_obs  = 200;
        sense_range = 6.0;
        w_vbar      = 200;
        if dist < d_safe_eff
            weight_obs = 2000;
        end
    end

    % v>=0 barrier
    if v < v_min + 0.5
        denom = (v - v_min + epsv);
        lx(3)    = lx(3) - (w_vbar / (denom^2));
        lxx(3,3) = lxx(3,3) + (2*w_vbar / (denom^3));
    end

    % obstacle barrier
    if dist < d_safe_eff + sense_range
        b = weight_obs * exp(-(dist - d_safe_eff));
        grad_dist = r / dist;
        lx(1:2) = lx(1:2) - b * grad_dist;

        H_unit = (eye(2)/dist) - (r*r')/(dist^3);
        lxx(1:2,1:2) = lxx(1:2,1:2) + b*(grad_dist*grad_dist' - H_unit);
        lxx(1:2,1:2) = 0.5*(lxx(1:2,1:2) + lxx(1:2,1:2)');
    end
end

function [lx, lxx] = terminal_cost(x, x_goal, y_target, v_target, Qf)
    dxT = [x(1)-x_goal; x(2)-y_target; x(3)-v_target; wrapToPiLocal(x(4))];
    lx  = Qf * dxT;
    lxx = Qf;
end

function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end

function J = compute_total_cost(x_traj, u_nom, obs_traj, Sigma_obs, ...
                                x_goal, y_target, v_target, ...
                                d_safe, beta, margin, stage, ...
                                Q, R, Qf, w_progress)
    N = size(u_nom,2);
    J = 0;

    for k = 1:N
        x = x_traj(:,k);
        u = u_nom(:,k);

        dx = [x(1)-x_goal; x(2)-y_target; x(3)-v_target; wrapToPiLocal(x(4))];

        % basic quadratic
        J = J + dx' * Q * dx + u' * R * u;

        % progress reward (linear shaping)
        J = J - w_progress * x(1);

        % inflated distance
        sigma_max  = sqrt(max(eig(Sigma_obs(:,:,k))));
        d_safe_eff = d_safe + beta*sigma_max + margin;

        % velocity barrier cost (match derivatives form)
        v = x(3);
        v_min = 0.05; epsv = 1e-3;
        if stage == 1
            w_vbar = 30;
            weight_obs = 30; sense_range = 4.0;
        else
            w_vbar = 200;
            weight_obs = 200; sense_range = 6.0;
        end
        if v < v_min + 0.5
            denom = (v - v_min + epsv);
            J = J + w_vbar / denom; % barrier (positive when close)
        end

        % obstacle barrier cost
        obs = obs_traj(:,k);
        dist = norm(x(1:2) - obs);
        if dist < d_safe_eff + sense_range
            b = weight_obs * exp(-(dist - d_safe_eff));
            J = J + b;
        end
    end

    % terminal
    xT = x_traj(:,end);
    dxT = [xT(1)-x_goal; xT(2)-y_target; xT(3)-v_target; wrapToPiLocal(xT(4))];
    J = J + dxT' * Qf * dxT;
end
