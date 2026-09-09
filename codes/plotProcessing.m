function fig = plotProcessing(filepath,showFigure)
%PLOTPROCESSING 无需阈值，显示真实 CSI 的处理阶段与逐窗曲线。
% addpath('codes'); plotProcessing('data/某个文件.csi');
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'codes','my_function'));
if nargin<1 || isempty(filepath), filepath=fullfile(root,'data','rx_2_260904_154612.csi'); end
if nargin<2, showFigure=true; end
C=configLoad; setupParser(C); [G,~]=fileLoad(filepath,C);
assert(~isempty(G),'RapidPD:NoCSI','没有可用的 HE-SU CSI。');
[~,gi]=max(arrayfun(@(g)numel(g.time_seconds),G)); g=G(gi);
w=winSplit(g.time_seconds,C); valid=find([w.valid]);
assert(~isempty(valid),'RapidPD:NoWindows','没有满足质量要求的完整窗口。');
F=numel(g.subcarrier_index); P=numel(g.time_seconds);
A=abs(double(reshape(g.csi_raw(1,1,:,:),P,F))); H=csiComplexNorm(A);
first=w(valid(1)); D=H(first.indices,:)-mean(H(first.indices,:),1);
mid=ceil(size(D,1)/2); [~,layers]=multiLayerACF(D(mid,:),C.acf_layers);
layerScores=nan(numel(w),C.acf_layers); score=nan(numel(w),1);
for wi=valid
    ids=w(wi).indices; n=numel(ids); acc=zeros(1,C.acf_layers);
    for rx=1:g.nrx
        for tx=1:g.ntx
            h=csiComplexNorm(reshape(g.csi_raw(tx,rx,ids,:),n,F));
            d=h-mean(h,1);
            if max(abs(d(:)))<=32*eps(max(h(:))), d(:)=0; end
            [~,ls]=multiLayerACF(d,C.acf_layers);
            for k=1:C.acf_layers, acc(k)=acc(k)+mean(ls{k}(:,2)); end
        end
    end
    layerScores(wi,:)=acc/(g.ntx*g.nrx); score(wi)=layerScores(wi,end);
end
visibility='off'; if showFigure, visibility='on'; end
fig=figure('Name','RapidPD - processing curves','Visible',visibility,'Color','w', ...
    'Position',[100 50 1300 930]); tiledlayout(3,2,'TileSpacing','compact');
tones=unique(round(linspace(1,F,5))); names=compose('sc %d',g.subcarrier_index(tones));
nexttile; plot(g.time_seconds,A(:,tones)); grid on;
title('1. Raw amplitude (Tx1-Rx1, five subcarriers)'); ylabel('|CSI|'); xlabel('Time (s)');
legend(names,'Location','best','NumColumns',2);
nexttile; plot(g.time_seconds,H(:,tones)); grid on;
title('2. Per-packet amplitude normalization'); ylabel('|CSI| / sum(|CSI|)'); xlabel('Time (s)');
nexttile; imagesc(g.subcarrier_index, g.time_seconds(first.indices), D); colorbar; axis xy;
title('3. Background-subtracted residual (first valid window)'); xlabel('Subcarrier index'); ylabel('Time (s)');
nexttile; hold on;
for k=1:C.acf_layers, plot(0:F-1,layers{k},'LineWidth',1); end
grid on; title('4. Subcarrier ACF (middle packet of first valid window)');
xlabel('Lag (retained subcarrier index)'); ylabel('Normalized ACF'); legend('Layer 1','Layer 2','Layer 3');
nexttile; plot([w.stop],layerScores,'LineWidth',1); grid on;
title('5. Lag-1 statistic at each layer (all streams averaged)'); xlabel('Time (s)'); ylabel('Window statistic');
legend('Layer 1','Layer 2','Layer 3','Location','best');
nexttile; plot([w.stop],score,'Color',[0.05 0.4 0.8],'LineWidth',1.2); grid on;
title('6. Final motion score (no threshold or class decision)'); xlabel('Time (s)'); ylabel('Motion statistic');
[~,name]=fileparts(filepath); sgtitle(string(name),'Interpreter','none');
out=fullfile(root,'results','processing'); if ~isfolder(out), mkdir(out); end
exportgraphics(fig,fullfile(out,[char(name) '_processing.png']),'Resolution',130);
T=table([w.stop].',score,layerScores(:,1),layerScores(:,2),layerScores(:,3), ...
    'VariableNames',{'time_s','score','layer1','layer2','layer3'});
writetable(T,fullfile(out,[char(name) '_curves.csv']));
fprintf('处理曲线：%s\n',fullfile(out,[char(name) '_processing.png']));
if ~showFigure, close(fig); end
end
