function fig = plotProcessing(filepath,showFigure)
%PLOTPROCESSING 显示 Desay FACF 各处理阶段及逐窗最终运动分数。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'codes','my_function'));
if nargin<1 || isempty(filepath), error('RapidPD:Input','请传入一份 .csi 文件路径。'); end
if nargin<2, showFigure=true; end
C=configLoad; setupParser(C); [G,~]=fileLoad(filepath,C);
assert(~isempty(G),'RapidPD:NoCSI','没有可用的 HE-SU CSI。');
[~,gi]=max(arrayfun(@(g)numel(g.time_seconds),G)); g=G(gi);
w=winSplit(g.time_seconds,C);
assert(~isempty(w),'RapidPD:NoWindows','没有完整的 20 包窗口。');

F=numel(g.subcarrier_index); P=numel(g.time_seconds);
A=abs(double(reshape(g.csi_raw(1,1,:,:),P,F)));
H=csiComplexNorm(A); first=w(1); h0=H(first.indices,:);
HH=[h0;h0]; dmin=C.pair_min_distance;
D=h0(1,:)-HH(1+dmin,:);
[~,layers]=desayFACF(h0,dmin,C.acf_shift);

score=nan(numel(w),1);
for wi=1:numel(w)
    d=winProcess(g.csi_raw(:,:,w(wi).indices,:),C);
    score(wi)=d.score;
end

visibility='off'; if showFigure, visibility='on'; end
fig=figure('Name','Desay FACF processing curves','Visible',visibility,'Color','w', ...
    'Position',[100 50 1300 930]); tiledlayout(3,2,'TileSpacing','compact');
tones=unique(round(linspace(1,F,5))); names=compose('sc %d',g.subcarrier_index(tones));
nexttile; plot(g.time_seconds,A(:,tones)); grid on;
title('1. Raw amplitude'); ylabel('|CSI|'); xlabel('Time (s)');
legend(names,'Location','best','NumColumns',2);
nexttile; plot(g.time_seconds,H(:,tones)); grid on;
title('2. Paper amplitude compensation'); ylabel('|CSI| / sum(|CSI|)'); xlabel('Time (s)');
nexttile; plot(g.subcarrier_index,D,'LineWidth',1); grid on;
title(sprintf('3. Packet-pair difference (packet 1 - packet %d)',1+dmin));
xlabel('Subcarrier index'); ylabel('Difference');
nexttile; imagesc(1:F-1,1:size(layers.first,1),layers.first); axis xy; colorbar;
title('4. First one-sided unbiased ACF'); xlabel('Positive lag'); ylabel('Packet pair');
nexttile; imagesc(1:F-2,1:size(layers.second,1),layers.second); axis xy; colorbar;
title('5. Second one-sided unbiased ACF'); xlabel('Positive lag'); ylabel('Packet pair');
nexttile; plot([w.stop],score,'Color',[0.05 0.4 0.8],'LineWidth',1.2); grid on;
title('6. Global normalized circular correlation'); xlabel('Time (s)'); ylabel('Motion statistic');
[~,name]=fileparts(filepath); sgtitle(string(name),'Interpreter','none');
out=fullfile(root,'results','processing'); if ~isfolder(out), mkdir(out); end
exportgraphics(fig,fullfile(out,[char(name) '_processing.png']),'Resolution',130);
T=table([w.start].',[w.stop].',repmat(C.window_packets,numel(w),1),score, ...
    'VariableNames',{'start_s','stop_s','packets','score'});
writetable(T,fullfile(out,[char(name) '_curves.csv']));
fprintf('处理曲线：%s\n',fullfile(out,[char(name) '_processing.png']));
if ~showFigure, close(fig); end
end
