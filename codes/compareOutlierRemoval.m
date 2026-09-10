function [Result,Summary,fig] = compareOutlierRemoval(dataA,dataB,varargin)
%COMPAREOUTLIERREMOVAL Compare FACF scores before/after whole-packet removal.
% A post-filter window must contain 20 uninterrupted original packets. It
% never bridges a removed packet or compresses the recording time axis.
assert(nargin>=2,'RapidPD:Inputs','Provide two parsed variables or two .csi paths.');
codeDir=fileparts(mfilename('fullpath')); addpath(fullfile(codeDir,'my_function'));
C=configLoad(varargin{:}); C.threshold=NaN;
inputs={dataA,dataB}; names=string({inputname(1),inputname(2)});
shifts=[1 4]; states=["before" "after"];
Result=struct('name',{},'filter',{},'before',{},'after',{},'metadata',{});
summaryRows={};
for j=1:2
    input=inputs{j};
    if ischar(input)||(isstring(input)&&isscalar(input))
        setupParser(C); [~,name]=fileparts(input); names(j)=string(name);
    elseif strlength(names(j))==0
        names(j)="Data "+j;
    end
    [groups,~]=fileLoad(input,C);
    assert(numel(groups)==1,'RapidPD:Groups','%s must contain one physical group.',names(j));
    g=groups(1); selectedRx=selectRxIndices(g.nrx,C);
    assert(numel(g.subcarrier_index)==1001 && all(diff(g.subcarrier_index)==1), ...
        'RapidPD:FrequencyGrid',[char(names(j)) ' must use the PicoScenes 1001-tone ' ...
        'uniform grid. Parse it with read_rxs_log before calling this experiment.']);
    csi=g.csi_raw(:,selectedRx,:,:);
    filter=packetSpectralOutliers(csi);
    fprintf('%s: removed %d/%d packets (%.2f%%), threshold %.6f.\n', ...
        names(j),sum(filter.removed),numel(filter.removed),100*mean(filter.removed),filter.threshold);
    for stateIndex=1:2
        if states(stateIndex)=="before"
            windows=winSplit(g.time_seconds,C);
        else
            windows=winSplitClean(g.time_seconds,filter.keep,C);
        end
        W=scoreWindows(csi,windows,C,shifts);
        Result(j).(states(stateIndex))=W;
        for si=1:numel(shifts)
            values=W.(sprintf('shift_%d',shifts(si)));
            summaryRows(end+1,:)={names(j),states(stateIndex),shifts(si), ... %#ok<AGROW>
                height(W),safeStat(values,@mean),safeStat(values,@std), ...
                safeStat(values,@median),safeStat(W.duration,@mean),safeStat(W.max_gap,@max)};
        end
    end
    metadata=rmfield(g,{'csi_raw','packet_csi'});
    metadata.selected_rx_indices=selectedRx;
    Result(j).name=names(j); Result(j).filter=filter; Result(j).metadata=metadata;
end
Summary=cell2table(summaryRows,'VariableNames', ...
    {'input','state','shift_bins','windows','mean_score','std_score', ...
    'median_score','mean_window_duration_s','max_interpacket_gap_s'});
disp(Summary);

visibility='off'; if C.show_figures, visibility='on'; end
fig=figure('Name','FACF before/after packet removal','Visible',visibility, ...
    'Color','w','Position',[80 80 1400 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for j=1:2
    for si=1:numel(shifts)
        nexttile; hold on;
        for stateIndex=1:2
            W=Result(j).(states(stateIndex));
            style='-'; if stateIndex==2, style='--'; end
            plot(W.stop,W.(sprintf('shift_%d',shifts(si))),style,'LineWidth',1.2, ...
                'DisplayName',states(stateIndex));
        end
        grid on; xlabel('Time from recording start (s)'); ylabel('Motion statistics');
        title(sprintf('%s, shift=%d',names(j),shifts(si)),'Interpreter','none');
        legend('Location','best');
    end
end
out=fullfile(C.path_res,['outlier_removal_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
mkdir(out);
for j=1:2
    writetable(Result(j).filter.table,fullfile(out,sprintf('input_%d_packet_filter.csv',j)));
    writetable(Result(j).before,fullfile(out,sprintf('input_%d_before.csv',j)));
    writetable(Result(j).after,fullfile(out,sprintf('input_%d_after.csv',j)));
end
writetable(Summary,fullfile(out,'summary.csv'));
save(fullfile(out,'outlier_removal_results.mat'),'Result','Summary','C','shifts','-v7.3');
exportgraphics(fig,fullfile(out,'motion_scores_before_after.png'),'Resolution',170);
fprintf('Output: %s\n',out);
if ~C.show_figures, close(fig); end
end

function W=scoreWindows(csi,timeWindows,C,shifts)
N=numel(timeWindows); start=nan(N,1); stop=start; duration=start; max_gap=start;
score=nan(N,numel(shifts));
for wi=1:N
    idx=timeWindows(wi).indices;
    start(wi)=timeWindows(wi).start; stop(wi)=timeWindows(wi).stop;
    duration(wi)=stop(wi)-start(wi); max_gap(wi)=timeWindows(wi).max_gap;
    for si=1:numel(shifts)
        Cs=C; Cs.acf_shift=shifts(si);
        d=winProcess(csi(:,:,idx,:),Cs); score(wi,si)=d.score;
    end
end
W=table(start,stop,duration,max_gap,repmat(C.window_packets,N,1), ...
    'VariableNames',{'start','stop','duration','max_gap','packets'});
for si=1:numel(shifts), W.(sprintf('shift_%d',shifts(si)))=score(:,si); end
end

function value=safeStat(x,fun)
x=x(~isnan(x));
if isempty(x), value=NaN; else, value=fun(x); end
end
