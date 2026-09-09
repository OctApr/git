function fig = plotComparison(resultsFile,labelsFile,profile,showFigure)
%PLOTCOMPARISON 显示与论文 Fig.6 同型的分数对比；没有标签时只用文件名。
% plotComparison() 自动打开最近一次结果，不必重跑原始数据。
% plotComparison(resultsMat,'data/labels.csv') 生成有真实标签的 Static/Motion 图。
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1 || isempty(resultsFile)
    files=dir(fullfile(root,'results','*','results.mat'));
    assert(~isempty(files),'RapidPD:NoResults','请先运行 main。');
    [~,i]=max([files.datenum]); resultsFile=fullfile(files(i).folder,files(i).name);
end
if nargin<2, labelsFile=''; end
if nargin<3, profile=''; end
if nargin<4, showFigure=true; end
S=load(resultsFile,'Results'); curves=struct('name',{},'profile',{},'score',{});
% 用户指定的默认对比文件；不按源 MAC 排除其中任意一份。
if isempty(labelsFile) && isempty(profile)
    names=strings(numel(S.Results),1);
    for fi=1:numel(S.Results)
        [~,name,ext]=fileparts(S.Results{fi}.file); names(fi)=string(name)+string(ext);
    end
    pair=["rx_2_260904_171616.csi","rx_2_260904_180803.csi"];
    if all(ismember(pair,names))
        fig=compareRecordings(resultsFile,pair,showFigure); return;
    end
end
for fi=1:numel(S.Results)
    R=S.Results{fi};
    for gi=1:numel(R.groups)
        W=R.groups{gi}.windows; ok=W.valid & isfinite(W.score);
        if ~any(ok), continue; end
        curves(end+1)=struct('name',string(R.file),'profile',string(W.profile(find(ok,1))), ...
            'score',W.score(ok)); %#ok<AGROW>
    end
end
assert(~isempty(curves),'RapidPD:NoResults','没有有效分数。');
if isempty(profile)
    p=string({curves.profile}); u=unique(p); counts=arrayfun(@(i)sum(p==u(i)),1:numel(u));
    [~,i]=max(counts); profile=u(i);
end
curves=curves(string({curves.profile})==string(profile));
assert(~isempty(curves),'RapidPD:Profile','找不到此 profile 的结果。');
visibility='off'; if showFigure, visibility='on'; end
fig=figure('Name','RapidPD - score comparison','Visible',visibility,'Color','w', ...
    'Position',[100 100 1150 720]);
tiledlayout(2,1,'TileSpacing','compact'); nexttile;
if ~isempty(labelsFile)
    T=buildLabeledWindows(resultsFile,labelsFile);
    T=T(T.profile==string(profile),:);
    if ~any(T.label==0) || ~any(T.label==1)
        close(fig); error('RapidPD:Labels','此 profile 缺少 Static/Motion 两类标签，不能绘制两类比较。');
    end
    plot(T.score(T.label==0),'-b','LineWidth',1.1); hold on;
    plot(T.score(T.label==1),'--','Color',[0.9 0.25 0.05],'LineWidth',1.1);
    legend('Static (label=0)','Motion/presence (label=1)','Location','best');
    title('Motion Statistics in Subcarrier Dimension - labeled recordings');
else
    % 同一 profile 的首尾两份文件；不根据分数猜测真实状态。
    chosen=unique([1 numel(curves)],'stable'); names=strings(1,numel(chosen));
    styles={'-','--'};
    for j=1:numel(chosen)
        c=curves(chosen(j)); plot(c.score,styles{j},'LineWidth',1.1); hold on;
        [~,name]=fileparts(c.name); names(j)=name;
    end
    legend(names,'Interpreter','none','Location','best');
    title('Motion Statistics - recording comparison (ground-truth labels unavailable)');
end
xlabel('Window index (valid windows)'); ylabel('Motion statistics'); grid on;
nexttile; mu=arrayfun(@(c)mean(c.score),curves);
lo=arrayfun(@(c)min(c.score),curves); hi=arrayfun(@(c)max(c.score),curves);
errorbar(1:numel(curves),mu,mu-lo,hi-mu,'o','LineWidth',1); grid on;
xlabel('Recording index (same source and radio profile)'); ylabel('Mean and min/max score');
title(string(profile),'Interpreter','none');
outputDir=fileparts(resultsFile);
exportgraphics(fig,fullfile(outputDir,'comparison.png'),'Resolution',140);
names=string({curves.name}).'; recording_index=(1:numel(curves)).';
writetable(table(recording_index,names,mu(:),lo(:),hi(:), ...
    'VariableNames',{'recording_index','file','mean_score','min_score','max_score'}), ...
    fullfile(outputDir,'comparison_recordings.csv'));
fprintf('对比图已生成：%s\n',fullfile(outputDir,'comparison.png'));
if ~showFigure, close(fig); end
end
