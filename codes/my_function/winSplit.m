function Windows = winSplit(t,Config)
%WINSPLIT 按 Desay 方式，以固定数量的连续数据包分窗。
% 时间戳仅用于绘图和记录；尾部不足 window_packets 的数据不进入计算。
t=t(:); Windows=struct('start',{},'stop',{},'indices',{},'max_gap',{},'valid',{});
if numel(t)<Config.window_packets, return; end
assert(all(diff(t)>0),'RapidPD:Time','时间戳必须严格递增。');
begins=1:Config.step_packets:(numel(t)-Config.window_packets+1);
for k=1:numel(begins)
    selected=begins(k)+(0:Config.window_packets-1);
    gaps=diff(t(selected));
    Windows(k)=struct('start',t(selected(1)),'stop',t(selected(end)), ...
        'indices',selected,'max_gap',max(gaps),'valid',true);
end
end
