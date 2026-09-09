function Windows = winSplit(t,Config)
%WINSPLIT 按真实时间分窗；对以首包为起点的 1/fs 网格选最近原始包。
% 不插值、不复制缺包；尾部不足完整窗口不生成阴性判决。
t=t(:); Windows=struct('start',{},'stop',{},'indices',{},'max_gap',{},'valid',{});
if numel(t)<2, return; end
assert(all(diff(t)>0),'RapidPD:Time','时间戳必须严格递增。');
dt=min(median(diff(t)),1/Config.target_rate);
captureEnd=t(end)+dt;
starts=t(1):Config.step_seconds:(captureEnd-Config.window_seconds+1e-6);
for k=1:numel(starts)
    a=starts(k); b=a+Config.window_seconds;
    candidates=find(t>=a & t<b);
    selected=[];
    if ~isempty(candidates)
        slot=floor((t(candidates)-a)*Config.target_rate+0.5);
        keep=slot<floor(Config.window_seconds*Config.target_rate);
        slot=slot(keep); candidates=candidates(keep);
        slots=unique(slot);
        for j=1:numel(slots)
            ids=candidates(slot==slots(j));
            center=a+slots(j)/Config.target_rate;
            [~,best]=min(abs(t(ids)-center)); selected(end+1)=ids(best); %#ok<AGROW>
        end
    end
    gaps=diff([a;t(selected);b]); maxGap=max(gaps);
    Windows(k)=struct('start',a,'stop',b,'indices',selected, ...
        'max_gap',maxGap,'valid',numel(selected)>=Config.min_packets && maxGap<=Config.max_gap_seconds);
end
end
