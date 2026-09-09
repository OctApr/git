function Summary = auditRapidPD
%AUDITRAPIDPD 固定数据的诊断实验，不修改检测算法或选择分类阈值。
% 用户确认：171616=一人运动；180258/180803/180942=无人。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'codes','my_function'));
C=configLoad; setupParser(C); rng(20260908);
files=["csi (1)一人/rx_2_260904_171616.csi", ...
    "csi(2)二人+无人/rx_2_260904_180258.csi", ...
    "csi(2)二人+无人/rx_2_260904_180803.csi", ...
    "csi(2)二人+无人/rx_2_260904_180942.csi"];
rows=cell(0,11); windowTables={};
for fi=1:numel(files)
    for interpolate=[false true]
        path=fullfile(root,'data',files(fi));
        if interpolate
            fid=fopen(path,'rb','ieee-le'); cleanup=onCleanup(@()fclose(fid));
            cells={};
            while true
                head=fread(fid,4,'*uint8'); if isempty(head), break; end
                n=double(typecast(head,'uint32')); body=fread(fid,n,'*uint8');
                cells{end+1}=RXSParser([head;body],true); %#ok<AGROW>
            end
            clear cleanup;
            [groups,audit]=fileLoad(cells,C);
        else
            [groups,audit]=fileLoad(path,C);
        end
        assert(numel(groups)==1 && audit.corrupt==0);
        g=groups(1); wins=winSplit(g.time_seconds,C); F=numel(g.subcarrier_index);
        shuffle=randperm(F); stats=nan(numel(wins),9);
        for wi=1:numel(wins)
            if ~wins(wi).valid, continue; end
            ids=wins(wi).indices; accum=zeros(1,9);
            for rx=1:g.nrx
                for tx=1:g.ntx
                    h=csiComplexNorm(reshape(g.csi_raw(tx,rx,ids,:),numel(ids),F));
                    d=h-mean(h,1); [~,ls]=multiLayerACF(d,3);
                    [~,sh]=multiLayerACF(d(:,shuffle),3);
                    h4=csiComplexNorm(reshape(g.csi_raw(tx,rx,ids,1:4:end),numel(ids),[]));
                    [~,s4]=multiLayerACF(h4-mean(h4,1),3);
                    for k=1:3, accum(k)=accum(k)+mean(ls{k}(:,2)); end
                    accum(4)=accum(4)+mean(sh{1}(:,2));
                    accum(5)=accum(5)+mean(sh{3}(:,2));
                    accum(6)=accum(6)+mean(s4{1}(:,2));
                    accum(7)=accum(7)+mean(s4{3}(:,2));
                    accum(8)=accum(8)+sqrt(mean(d(:).^2));
                    accum(9)=accum(9)+F*sqrt(mean(d(:).^2));
                end
            end
            stats(wi,:)=accum/(g.ntx*g.nrx);
        end
        m=mean(stats,1,'omitnan');
        rows(end+1,:)={files(fi),interpolate,F,sum(isfinite(stats(:,1))), ...
            m(1),m(2),m(3),m(4),m(5),m(6),m(7)}; %#ok<AGROW>
        T=array2table(stats,'VariableNames',{'layer1','layer2','layer3','shuffle_layer1', ...
            'shuffle_layer3','stride4_layer1','stride4_layer3','residual_rms','relative_residual_rms'});
        T.stop=[wins.stop].'; T.file=repmat(files(fi),height(T),1); T.interpolated=repmat(interpolate,height(T),1);
        windowTables{end+1}=T; %#ok<AGROW>
        fprintf('%s interp=%d L1=%.6f L2=%.6f L3=%.6f shuffleL3=%.6f\n',files(fi),interpolate,m(1:3),m(5));
    end
end
Summary=cell2table(rows,'VariableNames',{'file','interpolated','tones','windows', ...
    'layer1','layer2','layer3','shuffle_layer1','shuffle_layer3','stride4_layer1','stride4_layer3'});
out=fullfile(root,'results','audit'); if ~isfolder(out), mkdir(out); end
writetable(Summary,fullfile(out,'summary.csv')); Windows=vertcat(windowTables{:});
writetable(Windows,fullfile(out,'windows.csv'));
% Independent white-noise control: does three-layer ACF inevitably return ~1?
control=zeros(100,3);
for i=1:100
    h=csiComplexNorm(10+randn(20,980));
    [~,ls]=multiLayerACF(h-mean(h,1),3);
    for k=1:3, control(i,k)=mean(ls{k}(:,2)); end
end
NoiseControl=array2table(control,'VariableNames',{'layer1','layer2','layer3'});
writetable(NoiseControl,fullfile(out,'white_noise_control.csv'));
save(fullfile(out,'audit.mat'),'Summary','Windows','NoiseControl','C');
fprintf('WHITE NOISE mean L1=%.6f L2=%.6f L3=%.6f\n',mean(control));
f=figure('Visible','off','Color','w','Position',[100 100 1100 700]);
tiledlayout(2,1);
nexttile; values=[Summary.layer1 Summary.layer2 Summary.layer3];
bar(values(1:2:end,:)); grid on; ylim([0 1]);
xticklabels({'171616 motion','180258 empty','180803 empty','180942 empty'});
ylabel('Mean lag-1 score'); title('Actual files: layers 1, 2, 3 (interpolation OFF)');
legend('Layer 1','Layer 2','Layer 3','Location','best');
nexttile; bar([Summary.layer3(1:2:end),Summary.layer3(2:2:end),Summary.shuffle_layer3(1:2:end)]);
grid on; ylim([0 1]); xticklabels({'171616 motion','180258 empty','180803 empty','180942 empty'});
ylabel('Mean layer-3 score'); title('Controlled comparisons - no threshold fitted');
legend('Interpolation OFF','Interpolation ON','Randomized tone order (control)','Location','best');
exportgraphics(f,fullfile(out,'diagnosis.png'),'Resolution',140); close(f);
end
