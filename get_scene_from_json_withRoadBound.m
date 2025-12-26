function scene = get_scene_from_json_withRoadBound(jsonPath, sampleIdx, road_bound_showFlag, doPlot)
%GET_SCENE_FROM_JSON_WITHROADBOUND
% 从 ex_list_sample_updated.json 读取一个样本，提取：
% - lane_centerlines: cell，每条是 [Ni x 2] 的 (x_rel, y_rel)
% - road_boundaries:  struct/cell（取决于你 JSON 存法），原样返回
% - nv_past_xy:       [Nh x 2] ego历史(含当前点) (x_rel, y_rel)
% - nv_future_xy:     [Nf x 2] ego未来GT labels (x_rel, y_rel)
% - 以及必要的原始字段
%
% road_bound_showFlag: 1/0, 控制绘图是否显示道路边界
% doPlot: 1/0, 是否绘图（用于检查）

if nargin < 3, road_bound_showFlag = true; end
if nargin < 4, doPlot = true; end

txt  = fileread(jsonPath);
data = jsondecode(txt);

% 兼容：外层可能是数组/struct
samples = data;
if iscell(samples)
    m = samples{sampleIdx};
else
    m = samples(sampleIdx);
end

% 取核心字段
M     = m.matrix;                 % [L x 32]
spans = m.polyline_spans;         % slice-like struct array: start/stop
k0    = m.map_start_polyline_idx; % 前 k0 条 polyline 为 agents(ego+others)，之后为 lanes

% spans -> MATLAB 索引（cell）
S = cell(numel(spans),1);
for j = 1:numel(spans)
    st0 = spans(j).start;  % 0-based
    ed0 = spans(j).stop;   % open interval
    st = st0 + 1;          % MATLAB 1-based
    ed = ed0;              % MATLAB st:ed 对应 Python [st0, ed0-1]
    if ed < st
        S{j} = [];
    else
        S{j} = st:ed;
    end
end

agent_spans = S(1:k0);
lane_spans  = S(k0+1:end);

% ego (NV) past：约定 agent_spans{1} 是 ego
ego_span = agent_spans{1};
nv_past_xy = M(ego_span, 1:2);

% ego future GT labels
nv_future_xy = [];
if isfield(m, 'labels') && ~isempty(m.labels)
    nv_future_xy = m.labels;    % [Nf x 2] rel
end

% lane centerlines：x_curr=v[-3]=col30, y_curr=v[-4]=col29
lane_centerlines = cell(numel(lane_spans),1);
for l = 1:numel(lane_spans)
    idx = lane_spans{l};
    if isempty(idx), continue; end
    x_curr = M(idx, 30);
    y_curr = M(idx, 29);
    lane_centerlines{l} = [x_curr, y_curr];
end

% road boundaries（你新增的变量，原样拿出来）
road_boundaries = [];
if isfield(m, 'road_boundaries')
    road_boundaries = m.road_boundaries;
end

% 打包输出
scene = struct();
scene.sample = m;
scene.matrix = M;
scene.k0 = k0;
scene.lane_centerlines = lane_centerlines;
scene.road_boundaries  = road_boundaries;
scene.nv_past_xy = nv_past_xy;
scene.nv_future_xy = nv_future_xy;

% 可选绘图（论文风格：灰线+黑past+绿future）
if doPlot
    figure('Name','Scene (NV/Ego only, Local frame)'); hold on; grid on; axis equal;

    % lanes
    for l = 1:numel(lane_centerlines)
        xy = lane_centerlines{l};
        if isempty(xy), continue; end
        plot(xy(:,1), xy(:,2), '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0);
    end

    % road boundaries (optional)
    if road_bound_showFlag && ~isempty(road_boundaries)
        plot_road_boundaries_local(road_boundaries);
    end

    % NV past
    hPast = plot(nv_past_xy(:,1), nv_past_xy(:,2), 'k-', 'LineWidth', 2.5);
    hNow  = plot(nv_past_xy(end,1), nv_past_xy(end,2), 'ko', 'MarkerSize', 7, 'LineWidth', 2);

    % NV future
    hFut = [];
    hEnd = [];
    if ~isempty(nv_future_xy)
        hFut = plot(nv_future_xy(:,1), nv_future_xy(:,2), 'g-', 'LineWidth', 2.5);
        hEnd = plot(nv_future_xy(end,1), nv_future_xy(end,2), 'gs', 'MarkerSize', 7, 'LineWidth', 2);
    end

    xlabel('x_{rel} (m)'); ylabel('y_{rel} (m)');
    title(sprintf('Sample #%d', sampleIdx));

    % legend 用 handle 控制顺序（更稳）
    if isempty(hFut)
        legend([hPast, hNow], {'NV past','NV current'}, 'Location','bestoutside');
    else
        legend([hPast, hNow, hFut, hEnd], {'NV past','NV current','NV future (GT)','NV future end'}, ...
            'Location','bestoutside');
    end
end
end

function plot_road_boundaries_local(road_boundaries)
% 这里尽量“兼容”不同 JSON 存法：
% - 如果你存的是 cell，每个元素是 [N x 2]
% - 或者 struct array，字段可能叫 x/y 或 points
try
    if iscell(road_boundaries)
        for i = 1:numel(road_boundaries)
            bd = road_boundaries{i};
            if isempty(bd), continue; end
            plot(bd(:,1), bd(:,2), '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.8);
        end
    elseif isstruct(road_boundaries)
        for i = 1:numel(road_boundaries)
            rb = road_boundaries(i);
            if isfield(rb,'points')
                P = rb.points;
                plot(P(:,1), P(:,2), '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.8);
            elseif isfield(rb,'x') && isfield(rb,'y')
                plot(rb.x, rb.y, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.8);
            end
        end
    else
        % 如果直接是 numeric [N x 2]
        if isnumeric(road_boundaries) && size(road_boundaries,2)==2
            plot(road_boundaries(:,1), road_boundaries(:,2), '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.8);
        end
    end
catch
    % 失败就不画（不影响主流程）
end
end
