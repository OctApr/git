function [Result,figCurve,figHeatmap] = plotRapidPDFirstLayerACF(dataA,dataB,varargin)
%PLOTRAPIDPDFIRSTLAYERACF Plot the complete first-layer ACF for two inputs.
% Processing before ACF matches runRapidPDPaper: remove spectral-shape
% outlier packets, repack retained packets into 20-packet windows, normalize
% amplitude, calculate the per-tone window mean, then subtract that mean.
assert(nargin>=2,'RapidPD:Inputs','Provide two parsed variables or two .csi paths.');
codeDir=fileparts(mfilename('fullpath')); addpath(fullfile(codeDir,'my_function'));
C=configLoad(varargin{:});
inputs={dataA,dataB}; names=string({inputname(1),inputname(2)});
Result=struct('name',{},'filter',{},'streams',{},'metadata',{});

for j=1:2
    input=inputs{j};
    if ischar(input)||(isstring(input)&&isscalar(input))
        setupParser(C); [~,name]=fileparts(input); names(j)=string(name);
    elseif strlength(names(j))==0
        names(j)="Data "+j;
    end
    [groups,audit]=fileLoad(input,C);
    assert(numel(groups)==1,'RapidPD:Groups','%s must contain one physical group.',names(j));
    g=groups(1); selectedRx=selectRxIndices(g.nrx,C);
    csi=g.csi_raw(:,selectedRx,:,:);
    filter=packetSpectralOutliers(csi);
    kept=find(filter.keep); windowCount=floor(numel(kept)/C.window_packets);
    streamCount=size(csi,1)*size(csi,2); F=size(csi,4);
    streams=struct('window_mean_acf',{},'mean_acf',{},'std_acf',{});
    for s=1:streamCount
        streams(s).window_mean_acf=nan(windowCount,F);
        streams(s).mean_acf=nan(1,F);
        streams(s).std_acf=nan(1,F);
    end
    allRows=cell(1,streamCount);
    for s=1:streamCount, allRows{s}=nan(windowCount*C.window_packets,F); end
    for wi=1:windowCount
        idx=kept((wi-1)*C.window_packets+(1:C.window_packets));
        d=rapidPDPaperWindow(csi(:,:,idx,:),'acf_layers',1,'acf_lag',C.acf_shift);
        rows=(wi-1)*C.window_packets+(1:C.window_packets);
        for s=1:streamCount
            rho1=d.stream_details{s}.acf_layers{1};
            allRows{s}(rows,:)=rho1;
            streams(s).window_mean_acf(wi,:)=mean(rho1,1);
        end
    end
    for s=1:streamCount
        streams(s).mean_acf=mean(allRows{s},1,'omitnan');
        streams(s).std_acf=std(allRows{s},0,1,'omitnan');
    end
    metadata=rmfield(g,{'csi_raw','packet_csi'});
    metadata.selected_rx_indices=selectedRx; metadata.audit=audit;
    metadata.window_count=windowCount; metadata.acf_layer=1;
    Result(j)=struct('name',names(j),'filter',filter,'streams',streams,'metadata',metadata);
    fprintf('%s: removed %d/%d packets; %d windows; %d ACF lags.\n', ...
        names(j),sum(filter.removed),numel(filter.removed),windowCount,F);
end

assert(numel(Result(1).streams)==numel(Result(2).streams), ...
    'RapidPD:Streams','The two inputs must have the same selected stream count.');
streamCount=numel(Result(1).streams); visibility='off';
if C.show_figures, visibility='on'; end
styles={'-','--'};

figCurve=figure('Name','RapidPD complete first-layer ACF','Visible',visibility, ...
    'Color','w','Position',[70 80 1400 540]);
tiledlayout(1,streamCount,'TileSpacing','compact','Padding','compact');
for s=1:streamCount
    nexttile; hold on;
    for j=1:2
        lag=0:numel(Result(j).streams(s).mean_acf)-1;
        plot(lag,Result(j).streams(s).mean_acf,styles{j},'LineWidth',1.35, ...
            'DisplayName',Result(j).name);
    end
    yline(0,':','HandleVisibility','off'); grid on;
    xlabel('Lag (subcarrier index)'); ylabel('First-layer normalized ACF \rho_1');
    title(sprintf('Rx stream %d: complete ACF',s));
    legend('Interpreter','none','Location','best');
end
sgtitle('RapidPD first-layer ACF after window-mean subtraction');

figHeatmap=figure('Name','RapidPD first-layer ACF by window','Visible',visibility, ...
    'Color','w','Position',[60 60 1450 820]);
tiledlayout(2,streamCount,'TileSpacing','compact','Padding','compact');
for j=1:2
    for s=1:streamCount
        nexttile; imagesc(0:size(Result(j).streams(s).window_mean_acf,2)-1, ...
            1:size(Result(j).streams(s).window_mean_acf,1), ...
            Result(j).streams(s).window_mean_acf);
        axis xy; colorbar; xlabel('Lag (subcarrier index)'); ylabel('20-packet window');
        title(sprintf('%s | Rx stream %d',Result(j).name,s),'Interpreter','none');
    end
end
sgtitle('Mean first-layer ACF of the 20 packets in each window');

out=fullfile(C.path_res,['rapidpd_first_layer_acf_' ...
    char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]); mkdir(out);
for j=1:2
    writetable(Result(j).filter.table,fullfile(out,sprintf('input_%d_packet_filter.csv',j)));
    for s=1:streamCount
        lag=(0:numel(Result(j).streams(s).mean_acf)-1).';
        mean_acf=Result(j).streams(s).mean_acf.';
        std_acf=Result(j).streams(s).std_acf.';
        writetable(table(lag,mean_acf,std_acf), ...
            fullfile(out,sprintf('input_%d_stream_%d_complete_acf.csv',j,s)));
        writematrix(Result(j).streams(s).window_mean_acf, ...
            fullfile(out,sprintf('input_%d_stream_%d_window_acf.csv',j,s)));
    end
end
save(fullfile(out,'first_layer_acf_results.mat'),'Result','C','-v7.3');
exportgraphics(figCurve,fullfile(out,'complete_first_layer_acf.png'),'Resolution',180);
exportgraphics(figHeatmap,fullfile(out,'first_layer_acf_by_window.png'),'Resolution',180);
fprintf('Output: %s\n',out);
if ~C.show_figures, close(figCurve); close(figHeatmap); end
end
