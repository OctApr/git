function fig = compareRecordings(resultsFile,selectedFiles,showFigure)
%COMPARERECORDINGS 按文件名选取两份已有运动分数，不施加阈值。
% compareRecordings() 默认比较 171616 和 180803。
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1 || isempty(resultsFile)
    files=dir(fullfile(root,'results','*','results.mat'));
    assert(~isempty(files),'RapidPD:NoResults','请先运行 main。');
    [~,i]=max([files.datenum]); resultsFile=fullfile(files(i).folder,files(i).name);
end
if nargin<2 || isempty(selectedFiles)
    selectedFiles=["rx_2_260904_171616.csi","rx_2_260904_180803.csi"];
end
if nargin<3, showFigure=true; end
selectedFiles=string(selectedFiles);
assert(numel(selectedFiles)==2,'RapidPD:Selection','请选择两个文件。');
S=load(resultsFile,'Results'); scores=cell(1,2); profiles=strings(2,1);
for j=1:2
    matches={};
    for fi=1:numel(S.Results)
        R=S.Results{fi}; [~,name,ext]=fileparts(R.file);
        if string(name)+string(ext)~=selectedFiles(j), continue; end
        for gi=1:numel(R.groups)
            W=R.groups{gi}.windows;
            if any(W.valid & isfinite(W.score)), matches{end+1}=W; end %#ok<AGROW>
        end
    end
    assert(numel(matches)==1,'RapidPD:Selection', ...
        '%s 应匹配一个有效分组，实际为 %d。',selectedFiles(j),numel(matches));
    W=matches{1}; ok=W.valid & isfinite(W.score);
    scores{j}=W.score(ok); profiles(j)=W.profile(find(ok,1));
end
visibility='off'; if showFigure, visibility='on'; end
fig=figure('Name','RapidPD - 171616 vs 180803','Visible',visibility, ...
    'Color','w','Position',[120 120 1100 500]);
plot(1:numel(scores{1}),scores{1},'-','Color',[0 0.447 0.741],'LineWidth',1.4); hold on;
plot(1:numel(scores{2}),scores{2},'--','Color',[0.85 0.325 0.098],'LineWidth',1.4);
legend(erase(selectedFiles,'.csi'),'Interpreter','none','Location','best');
title('Motion Statistics in Subcarrier Dimension');
xlabel('Window index (valid windows)'); ylabel('Motion statistics'); grid on;
out=fileparts(resultsFile);
stem=char(join(erase(selectedFiles,'.csi'),'_vs_'));
exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',160);
exportgraphics(fig,fullfile(out,'comparison.png'),'Resolution',160);
N=max(cellfun(@numel,scores)); values=nan(N,2);
for j=1:2, values(1:numel(scores{j}),j)=scores{j}; end
T=table((1:N).',values(:,1),values(:,2),'VariableNames', ...
    {'window_index',char(erase(selectedFiles(1),'.csi')),char(erase(selectedFiles(2),'.csi'))});
writetable(T,fullfile(out,[stem '.csv']));
M=table(selectedFiles(:),profiles,cellfun(@numel,scores).', ...
    cellfun(@mean,scores).','VariableNames',{'file','profile','valid_windows','mean_score'});
writetable(M,fullfile(out,'comparison_recordings.csv'));
disp(M(:,{'file','valid_windows','mean_score'}));
fprintf('曲线已保存：%s\n',fullfile(out,[stem '.png']));
if ~showFigure, close(fig); end
end
