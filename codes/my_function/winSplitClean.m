function Windows = winSplitClean(t,keep,Config)
%WINSPLITCLEAN Form windows only inside uninterrupted runs of kept packets.
% Returned indices address the original CSI array. No window crosses a
% removed packet, so packet removal cannot silently stretch the time axis.
t=t(:); keep=logical(keep(:));
assert(numel(t)==numel(keep),'RapidPD:CleanWindows','Time and keep mask sizes differ.');
assert(all(diff(t)>0),'RapidPD:Time','Timestamps must be strictly increasing.');
Windows=struct('start',{},'stop',{},'indices',{},'max_gap',{},'valid',{});
edges=diff([false;keep;false]);
runStart=find(edges==1); runStop=find(edges==-1)-1;
for r=1:numel(runStart)
    lastStart=runStop(r)-Config.window_packets+1;
    if lastStart<runStart(r), continue; end
    for first=runStart(r):Config.step_packets:lastStart
        selected=first+(0:Config.window_packets-1);
        gaps=diff(t(selected));
        Windows(end+1)=struct('start',t(selected(1)),'stop',t(selected(end)), ... %#ok<AGROW>
            'indices',selected,'max_gap',max(gaps),'valid',true);
    end
end
end
