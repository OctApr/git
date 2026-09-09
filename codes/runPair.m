function [Result,fig] = runPair(dataA,dataB,varargin)
%RUNPAIR 仅处理两个输入：工作区已解析的变量，或者两个 .csi 路径。
% [R,fig]=runPair(rx_2_260904_171616,rx_2_260904_180803);
% [R,fig]=runPair('D:/a.csi','D:/b.csi');
% 可加 'show_figures',false 等 configLoad 参数；默认只输出分数。
assert(nargin>=2,'RapidPD:Inputs','请输入两份已解析变量，或两个 .csi 文件路径。');
codeDir=fileparts(mfilename('fullpath')); addpath(fullfile(codeDir,'my_function'));
C=configLoad(varargin{:}); C.threshold=NaN;
fprintf('幅度补偿：paper_normalization（AX210 帧无 Desay 等价 AGC 总增益字段）。\n');
inputs={dataA,dataB}; names=string({inputname(1),inputname(2)});
Result=struct('name',{},'audit',{},'metadata',{},'windows',{},'config',{});
for j=1:2
    input=inputs{j};
    if ischar(input) || (isstring(input)&&isscalar(input))
        setupParser(C); [~,name]=fileparts(input); names(j)=string(name);
    elseif strlength(names(j))==0
        names(j)="Data "+j;
    end
    [groups,audit]=fileLoad(input,C);
    if audit.interpolated_frames>0
        fprintf('%s：%d 个包包含官方插值/预处理结果，使用工作区现有 CSI。\n',names(j),audit.interpolated_frames);
    end
    assert(numel(groups)==1,'RapidPD:Groups', ...
        '%s 得到 %d 个物理分组；请用 source_mac/packet_formats 筛选或传入单一分组。',names(j),numel(groups));
    g=groups(1); w=winSplit(g.time_seconds,C); N=numel(w);
    start=nan(N,1); stop=start; score=start; packets=zeros(N,1); valid=false(N,1);
    stream=nan(N,g.ntx*g.nrx);
    for wi=1:N
        start(wi)=w(wi).start; stop(wi)=w(wi).stop;
        packets(wi)=numel(w(wi).indices); valid(wi)=w(wi).valid;
        if ~valid(wi), continue; end
        d=winProcess(g.csi_raw(:,:,w(wi).indices,:),C);
        score(wi)=d.score; stream(wi,:)=d.stream_statistics;
    end
    assert(any(valid),'RapidPD:NoWindows','%s 没有有效时间窗口。',names(j));
    W=table(start,stop,packets,valid,score);
    for si=1:size(stream,2), W.(sprintf('stream_%d',si))=stream(:,si); end
    Result(j)=struct('name',names(j),'audit',audit, ...
        'metadata',rmfield(g,{'csi_raw','packet_csi'}),'windows',W,'config',C);
    fprintf('%s：%d 包，%d 个有效窗口，平均运动分数 %.6f\n', ...
        names(j),size(g.csi_raw,3),sum(valid),mean(score,'omitnan'));
end
visibility='off'; if C.show_figures, visibility='on'; end
fig=figure('Name','Desay FACF - two inputs','Visible',visibility,'Color','w','Position',[100 100 1100 500]);
styles={'-','--'};
for j=1:2
    W=Result(j).windows;
    plot(W.stop,W.score,styles{j},'LineWidth',1.4); hold on;
end
legend(names,'Interpreter','none','Location','best'); grid on;
title('Desay FACF Motion Statistics (20 packets/window)');
xlabel('Time from each recording start (s)'); ylabel('Motion statistics');
out=fullfile(C.path_res,['pair_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]); mkdir(out);
for j=1:2, writetable(Result(j).windows,fullfile(out,sprintf('input_%d.csv',j))); end
save(fullfile(out,'pair_results.mat'),'Result','C','-v7.3');
exportgraphics(fig,fullfile(out,'motion_scores.png'),'Resolution',150);
fprintf('双数据结果：%s\n',out);
if ~C.show_figures, close(fig); end
end
