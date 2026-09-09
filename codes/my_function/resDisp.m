function resDisp(W,G,Config,path)
%RESDISP 离线诊断图；关闭图窗，避免批量运行堆积 UI。
f=figure('Visible','off','Color','w','Position',[100 100 1200 650]);
cleanup=onCleanup(@()close(f)); %#ok<NASGU>
tiledlayout(2+double(any(isfinite(W.raw))),1,'TileSpacing','compact');
nexttile; plot(W.stop,W.score,'LineWidth',1); hold on;
if isfinite(Config.threshold), yline(Config.threshold,'--r','threshold'); end
ylabel('Motion score'); grid on;
title(sprintf('RapidPD | %s | %d MHz | %d tones | %d Tx x %d Rx', ...
    G.source_mac,G.bandwidth_mhz,numel(G.subcarrier_index),G.ntx,G.nrx),'Interpreter','none');
if any(isfinite(W.raw))
    nexttile;
    stairs(W.stop,W.raw,':','LineWidth',1); hold on;
    stairs(W.stop,W.smooth,'LineWidth',1.3); ylim([-0.1 1.1]); ylabel('Presence');
    legend('Raw','Causal majority','Location','best'); grid on;
    xlim([min(W.start),max(W.stop)]);
end
nexttile; plot(W.stop,W.packets,'LineWidth',1); hold on;
yline(Config.window_packets,':r'); ylabel('Packets/window'); xlabel('Time from file start (s)'); grid on;
exportgraphics(f,path,'Resolution',130);
end
