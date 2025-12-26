%% plot_nuscenes_ego_only_with_boundaries_from_json.m
% 只绘制：
%   - Lane centerlines (gray)
%   - Road boundaries (dark gray / black)   <-- 新增：road_boundaries
%   - Ego past trajectory (black)
%   - Ego future GT labels (green)
% 不绘制任何邻车/行人轨迹
%
% -------------------------------------------------------------

clear; clc; close all;
road_bound_showFlag = true;   % true: 显示道路边界 | false: 只显示车道中心线

%% ========== 1) 读取 JSON ==========
jsonPath = 'ex_list_sample_updated.json';   % 放到当前 MATLAB 工作目录下
txt  = fileread(jsonPath);
data = jsondecode(txt);

% 如果你的 JSON 外层结构不是“样本数组”，改这里
samples = data;

% 选择第几个样本
i = 1;
if iscell(samples)
    m = samples{i};
else
    m = samples(i);
end

%% ========== 2) 取核心字段 ==========
M     = m.matrix;                 % [L x 32]
spans = m.polyline_spans;         % struct array: start/stop(/step)
k0    = m.map_start_polyline_idx; % 前 k0 条 polyline 是轨迹(ego+agents)，之后是车道

if isempty(M) || isempty(spans)
    error('matrix 或 polyline_spans 为空：请确认导出 JSON 字段。');
end

%% ========== 3) spans -> MATLAB 索引（cell） ==========
S = cell(numel(spans),1);
for j = 1:numel(spans)
    st0 = spans(j).start;  % 0-based
    ed0 = spans(j).stop;   % slice stop, 开区间
    st = st0 + 1;          % MATLAB 1-based
    ed = ed0;              % MATLAB 用 st:ed 对应 [st0, ed0-1]
    if ed < st
        S{j} = [];
    else
        S{j} = st:ed;
    end
end

agent_spans = S(1:k0);
lane_spans  = S(k0+1:end);

%% ========== 4) Ego 历史轨迹 ==========
% 约定：第 1 条 agent polyline 是 ego（你的 encode_agents 就是这样）
ego_span = agent_spans{1};
ego_past_xy = M(ego_span, 1:2);   % agent polyline 的 (x,y) 在 v[0],v[1]

%% ========== 5) Ego 未来 GT（labels, 局部系） ==========
ego_future_xy = [];
if isfield(m, 'labels') && ~isempty(m.labels)
    ego_future_xy = m.labels;     % [12x2] (x_rel, y_rel)
end

%% ========== 6) Lane polylines（灰色）==========
% 重要：encode_lanes() 里 i_point==0 被跳过
% 所以每条 lane span 的第一行存的是 (prev -> curr) 的第一条边
% 若只画 curr，会丢掉起点，导致画出来和之前不一致。
%
% lane row:
%   x_curr = v[-3] -> MATLAB col 30
%   y_curr = v[-4] -> MATLAB col 29
%   x_prev = v[-1] -> MATLAB col 32
%   y_prev = v[-2] -> MATLAB col 31

lane_xy = cell(numel(lane_spans),1);
for l = 1:numel(lane_spans)
    idx = lane_spans{l};
    if isempty(idx), continue; end

    % curr 序列
    x_curr = M(idx, 30);
    y_curr = M(idx, 29);

    % 起点用第一行的 prev（把丢掉的 i_point=0 补回来）
    x0 = M(idx(1), 32);
    y0 = M(idx(1), 31);

    lane_xy{l} = [x0, y0; x_curr, y_curr];
end


%% ========== 7) Road boundaries（新增） ==========
% 你新增的 road_boundaries 建议结构：
%   m.road_boundaries = { [N1x2], [N2x2], ... }   或者 JSON 解码后为 cell 嵌套
road_bd = {};
if isfield(m, 'road_boundaries') && ~isempty(m.road_boundaries)
    road_bd = normalizePolylineList(m.road_boundaries);
end

%% ========== 8) 绘图（论文风格一些） ==========
fig = figure('Name','nuScenes (Ego only, Local frame)');
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');

% 统一一些排版
set(ax, 'FontName','Times New Roman', 'FontSize',12, 'LineWidth',1.0);
xlabel(ax, 'x_{rel} (m)');
ylabel(ax, 'y_{rel} (m)');

h_lane = gobjects(0);
h_bd   = gobjects(0);
h_past = gobjects(0);
h_curr = gobjects(0);
h_fut  = gobjects(0);
h_fend = gobjects(0);

% --- 车道中心线：灰色（虚线更像 centerline）
for l = 1:numel(lane_xy)
    if isempty(lane_xy{l}), continue; end
    if isempty(h_lane)
        h_lane = plot(ax, lane_xy{l}(:,1), lane_xy{l}(:,2), '-', ...
            'Color',[0.65 0.65 0.65], 'LineWidth',1.2, 'DisplayName','Lane centerlines');
    else
        plot(ax, lane_xy{l}(:,1), lane_xy{l}(:,2), '-', ...
            'Color',[0.65 0.65 0.65], 'LineWidth',1.2, 'HandleVisibility','off');
    end
end

% --- 道路边界：深灰/黑色（实线）
% --- 道路边界：仅在开关打开时绘制 ---
if road_bound_showFlag && ~isempty(road_bd)
    for b = 1:numel(road_bd)
        xy = road_bd{b};
        if isempty(xy) || size(xy,2)~=2, continue; end

        if isempty(h_bd)
            h_bd = plot(ax, xy(:,1), xy(:,2), '-', ...
                'Color',[0.15 0.15 0.15], 'LineWidth',1.8, ...
                'DisplayName','Road boundaries');
        else
            plot(ax, xy(:,1), xy(:,2), '-', ...
                'Color',[0.15 0.15 0.15], 'LineWidth',1.8, ...
                'HandleVisibility','off');
        end
    end
end

% --- Ego 历史：黑色
h_past = plot(ax, ego_past_xy(:,1), ego_past_xy(:,2), 'k-', ...
    'LineWidth',2.6, 'DisplayName','Ego past');

% 当前点（历史末端）
h_curr = plot(ax, ego_past_xy(end,1), ego_past_xy(end,2), 'ko', ...
    'MarkerSize',7, 'LineWidth',1.8, 'DisplayName','Ego current');

% --- Ego 未来：绿色（GT labels）
if ~isempty(ego_future_xy)
    h_fut = plot(ax, ego_future_xy(:,1), ego_future_xy(:,2), 'g-', ...
        'LineWidth',2.6, 'DisplayName','Ego future (GT)');
    h_fend = plot(ax, ego_future_xy(end,1), ego_future_xy(end,2), 'gs', ...
        'MarkerSize',7, 'LineWidth',1.8, 'DisplayName','Ego future end');
end

title(ax, sprintf('Sample #%d (Local frame)', i));

% legend：让 MATLAB 根据 DisplayName 自动生成（你说的“直接根据线段和点决定标注”）
legend(ax, 'show', 'Location','bestoutside');

%% ========== 9) 输出基本信息 ==========
fprintf('--- Sample #%d ---\n', i);
fprintf('matrix: %d x %d\n', size(M,1), size(M,2));
fprintf('polylines: total=%d, agents=%d, lanes=%d\n', numel(S), numel(agent_spans), numel(lane_spans));
fprintf('ego past points: %d\n', size(ego_past_xy,1));
if ~isempty(ego_future_xy)
    fprintf('ego future points (labels): %d\n', size(ego_future_xy,1));
end
fprintf('road boundaries polylines: %d\n', numel(road_bd));

%% ===================== Helper =====================
function out = normalizePolylineList(x)
% 将 jsondecode 后的 “多条 polyline” 统一成 cell{K}，每个元素是 [N x 2] double
    out = {};
    if isempty(x), return; end

    % 情况1：已经是数值矩阵 [N x 2]
    if isnumeric(x) && size(x,2)==2
        out = {double(x)};
        return;
    end

    % 情况2：cell，每个 cell 里可能是 [N x 2] 或者更深层 cell
    if iscell(x)
        tmp = {};
        for i = 1:numel(x)
            yi = x{i};
            if isempty(yi), continue; end
            if isnumeric(yi) && size(yi,2)==2
                tmp{end+1,1} = double(yi); %#ok<AGROW>
            elseif iscell(yi)
                % 递归拍平
                sub = normalizePolylineList(yi);
                for k = 1:numel(sub)
                    tmp{end+1,1} = sub{k}; %#ok<AGROW>
                end
            end
        end
        out = tmp;
        return;
    end

    % 其他结构（struct等）：尽量兼容常见形式：[{x:..,y:..}, ...] 这种就需要你告诉字段名
    warning('road_boundaries 的结构不是 numeric/cell，当前脚本未覆盖该格式。');
end
