function [Summary,Result,fig,outputDir] = compareShifts(dataA,dataB,varargin)
%COMPARESHIFTS 在同一1001-tone输入上对比Desay FACF的shift=1与shift=4。
% 输入应为PicoScenes默认插值后拖入工作区的变量；这样相邻列对应
% 78.125 kHz，shift=4对应RapidPD设备约312.5 kHz的频率间隔。
assert(nargin>=2,'RapidPD:Inputs','请输入两份已解析的PicoScenes工作区变量。');
codeDir=fileparts(mfilename('fullpath')); addpath(fullfile(codeDir,'my_function'));
C=configLoad(varargin{:}); C.threshold=NaN; C.stream_aggregation='mean';
shifts=[1 4]; inputs={dataA,dataB}; names=string({inputname(1),inputname(2)});
Result=repmat(struct('name',"",'audit',struct(),'metadata',struct(),'windows',table()),2,1);
summaryRows=cell(0,7);

for j=1:2
    assert(~(ischar(inputs{j}) || (isstring(inputs{j})&&isscalar(inputs{j}))), ...
        'RapidPD:InterpolatedInput', ...
        '本实验要求1001-tone插值数据：请先用read_rxs_log读取或拖入工作区，再传变量。');
    if strlength(names(j))==0, names(j)="Data "+j; end
    [groups,audit]=fileLoad(inputs{j},C);
    assert(isscalar(groups),'RapidPD:Groups','%s包含多个物理分组。',names(j));
    g=groups(1);
    assert(g.bandwidth_mhz==80 && g.packet_format==3, 'RapidPD:PHY', ...
        '%s必须是80 MHz HE数据。',names(j));
    assert(numel(g.subcarrier_index)==1001 && ...
        isequal(g.subcarrier_index(:).',-500:500), 'RapidPD:ToneGrid', ...
        '%s不是-500:500的1001-tone均匀频率栅格。',names(j));
    selectedRx=selectRxIndices(g.nrx,C); wins=winSplit(g.time_seconds,C); N=numel(wins);
    start=nan(N,1); stop=start; packets=zeros(N,1); valid=false(N,1);
    score=nan(N,numel(shifts)); stream=nan(N,g.ntx*numel(selectedRx),numel(shifts));
    for wi=1:N
        start(wi)=wins(wi).start; stop(wi)=wins(wi).stop;
        packets(wi)=numel(wins(wi).indices); valid(wi)=wins(wi).valid;
        if ~valid(wi), continue; end
        block=g.csi_raw(:,selectedRx,wins(wi).indices,:);
        for si=1:numel(shifts)
            Cs=C; Cs.acf_shift=shifts(si);
            detail=winProcess(block,Cs);
            score(wi,si)=detail.score;
            stream(wi,:,si)=detail.stream_statistics;
        end
    end
    assert(any(valid),'RapidPD:NoWindows','%s没有有效的20包窗口。',names(j));
    W=table(start,stop,packets,valid,score(:,1),score(:,2),score(:,2)-score(:,1), ...
        'VariableNames',{'start','stop','packets','valid','shift_1','shift_4','shift4_minus_shift1'});
    for si=1:numel(shifts)
        for streamIdx=1:size(stream,2)
            W.(sprintf('stream_%d_shift_%d',streamIdx,shifts(si)))=stream(:,streamIdx,si);
        end
        values=score(valid,si);
        summaryRows(end+1,:)={names(j),shifts(si),78.125*shifts(si), ...
            sum(valid),mean(values),std(values),median(values)}; %#ok<AGROW>
    end
    metadata=rmfield(g,{'csi_raw','packet_csi'}); metadata.selected_rx_indices=selectedRx;
    Result(j)=struct('name',names(j),'audit',audit,'metadata',metadata,'windows',W);
end

Summary=cell2table(summaryRows,'VariableNames', ...
    {'input','shift_bins','frequency_lag_khz','valid_windows','mean_score','std_score','median_score'});
visibility='off'; if C.show_figures, visibility='on'; end
fig=figure('Name','AX210 shift=1 vs shift=4','Visible',visibility, ...
    'Color','w','Position',[100 100 1100 760]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact'); styles={'-','--'};
for si=1:numel(shifts)
    nexttile; hold on
    for j=1:2
        plot(Result(j).windows.stop,Result(j).windows.(sprintf('shift_%d',shifts(si))), ...
            styles{j},'LineWidth',1.3,'DisplayName',names(j));
    end
    grid on; legend('Interpreter','none','Location','best');
    ylabel('Motion statistics'); xlabel('Time from recording start (s)');
    title(sprintf('Desay FACF: shift=%d (%g kHz)',shifts(si),78.125*shifts(si)));
end

outputDir=fullfile(C.path_res,['shift_ablation_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
mkdir(outputDir); writetable(Summary,fullfile(outputDir,'summary.csv'));
for j=1:2, writetable(Result(j).windows,fullfile(outputDir,sprintf('input_%d.csv',j))); end
save(fullfile(outputDir,'shift_ablation.mat'),'Summary','Result','C','shifts','-v7.3');
exportgraphics(fig,fullfile(outputDir,'shift_comparison.png'),'Resolution',160);
disp(Summary); fprintf('移位量对照结果：%s\n',outputDir);
if ~C.show_figures, close(fig); end
end
