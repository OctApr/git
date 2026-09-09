function Config = configLoad(varargin)
%CONFIGLOAD RapidPD 参数。时间单位秒，窗口独立于实际采样率。
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
Config = struct('path_data',fullfile(root,'data'), ...
    'path_res',fullfile(root,'results'), 'parser_path','', ...
    'window_seconds',1, 'step_seconds',1, 'target_rate',20, ...
    'min_packets',16, 'max_gap_seconds',0.20, ...
    'acf_layers',3, 'acf_mode','recursive', ...
    'time_aggregation','mean', 'stream_aggregation','mean', ...
    'threshold',NaN, 'threshold_profile','', 'smooth_windows',3, 'make_plots',true, ...
    'show_figures',true, ...
    'file_pattern','*.csi','packet_formats',3,'source_mac','');
assert(mod(numel(varargin),2)==0,'RapidPD:Config','参数必须成对出现。');
for k=1:2:numel(varargin)
    key=char(varargin{k});
    assert(isfield(Config,key),'RapidPD:Config','未知参数: %s',key);
    Config.(key)=varargin{k+1};
end
validateattributes(Config.window_seconds,{'numeric'},{'scalar','positive','finite'});
validateattributes(Config.step_seconds,{'numeric'},{'scalar','positive','finite'});
validateattributes(Config.target_rate,{'numeric'},{'scalar','integer','positive'});
validateattributes(Config.min_packets,{'numeric'},{'scalar','integer','>=',2});
validateattributes(Config.max_gap_seconds,{'numeric'},{'scalar','positive','finite'});
validateattributes(Config.acf_layers,{'numeric'},{'scalar','integer','positive'});
validateattributes(Config.smooth_windows,{'numeric'},{'scalar','integer','positive'});
assert(mod(Config.smooth_windows,2)==1,'RapidPD:Config','投票窗口数必须为奇数。');
assert(Config.min_packets<=floor(Config.window_seconds*Config.target_rate), ...
    'RapidPD:Config','min_packets 超过窗口采样容量。');
assert(isscalar(Config.threshold) && (isnan(Config.threshold)||isfinite(Config.threshold)), ...
    'RapidPD:Config','threshold 必须是有限标量或 NaN。');
assert(ismember(string(Config.acf_mode),["recursive","literal"]),'RapidPD:Config','acf_mode 无效。');
assert(ismember(string(Config.time_aggregation),["mean","sum"]),'RapidPD:Config','time_aggregation 无效。');
assert(ismember(string(Config.stream_aggregation),["mean","sum"]),'RapidPD:Config','stream_aggregation 无效。');
end
