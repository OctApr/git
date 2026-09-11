function [Result,Summary,fig] = runRapidPDPaper(dataA,dataB,varargin)
%RUNRAPIDPDPAPER Run the paper benchmark-removal path on two AX210 inputs.
% Whole-packet spectral-shape outliers are removed first. Retained packets
% are then kept in acquisition order and repacked into fixed 20-packet
% windows before any formal amplitude normalization is performed.
assert(nargin>=2,'RapidPD:Inputs','Provide two parsed variables or two .csi paths.');
codeDir=fileparts(mfilename('fullpath')); addpath(fullfile(codeDir,'my_function'));
C=configLoad(varargin{:}); C.threshold=NaN;
inputs={dataA,dataB}; names=string({inputname(1),inputname(2)});
Result=struct('name',{},'filter',{},'windows',{},'metadata',{});
rows={};
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

    % Diagnostic copy only: packetSpectralOutliers normalizes internally to
    % identify spectral shape, but formal algorithm normalization occurs in
    % rapidPDPaperWindow after rejected packets have been deleted.
    filter=packetSpectralOutliers(csi);
    kept=find(filter.keep); N=floor(numel(kept)/20);
    start=nan(N,1); stop=start; max_gap=start; phi_mean=start;
    Phi_text_average=start; Phi_formula_sum=start; Phi_literal_formula_sum=start;
    stream_phi_average=nan(N,size(csi,1)*size(csi,2));
    stream_phi_formula_sum=nan(N,size(csi,1)*size(csi,2));
    for wi=1:N
        idx=kept((wi-1)*20+(1:20));
        d=rapidPDPaperWindow(csi(:,:,idx,:), ...
            'acf_layers',C.paper_acf_layers,'acf_lag',C.acf_shift,'acf_mode','recursive');
        literal=rapidPDPaperWindow(csi(:,:,idx,:), ...
            'acf_layers',C.paper_acf_layers,'acf_lag',C.acf_shift,'acf_mode','literal');
        start(wi)=g.time_seconds(idx(1)); stop(wi)=g.time_seconds(idx(end));
        max_gap(wi)=max(diff(g.time_seconds(idx)));
        phi_mean(wi)=d.diagnostic_stream_mean;
        Phi_text_average(wi)=d.Phi_text_average;
        Phi_formula_sum(wi)=d.Phi_formula_sum;
        Phi_literal_formula_sum(wi)=literal.Phi_formula_sum;
        stream_phi_average(wi,:)=d.stream_phi_average;
        stream_phi_formula_sum(wi,:)=d.stream_phi_formula_sum;
    end
    W=table(start,stop,max_gap,repmat(20,N,1),phi_mean,Phi_text_average, ...
        Phi_formula_sum,Phi_literal_formula_sum, ...
        'VariableNames',{'start','stop','max_gap','packets','diagnostic_stream_mean', ...
        'Phi_text_average','Phi_formula_sum','Phi_literal_formula_sum'});
    for s=1:size(stream_phi_formula_sum,2)
        W.(sprintf('stream_%d_phi_formula_sum',s))=stream_phi_formula_sum(:,s);
        W.(sprintf('stream_%d_phi_time_average',s))=stream_phi_average(:,s);
    end
    metadata=rmfield(g,{'csi_raw','packet_csi'});
    metadata.selected_rx_indices=selectedRx; metadata.audit=audit;
    metadata.processing_order={'remove_spectral_outliers','amplitude_normalize', ...
        'per_tone_window_mean','subtract_window_mean', ...
        sprintf('%d_layer_acf',C.paper_acf_layers),'lag_statistic'};
    Result(j)=struct('name',names(j),'filter',filter,'windows',W,'metadata',metadata);
    for s=1:size(stream_phi_formula_sum,2)
        rows(end+1,:)={names(j),s,size(csi,3),sum(filter.removed),sum(filter.keep),N, ... %#ok<AGROW>
            safeMean(stream_phi_formula_sum(:,s)),safeMean(stream_phi_average(:,s))};
        fprintf('%s, Rx stream %d: equation-(19) phi %.6f; time average %.6f\n', ...
            names(j),s,safeMean(stream_phi_formula_sum(:,s)),safeMean(stream_phi_average(:,s)));
    end
end
Summary=cell2table(rows,'VariableNames',{'input','stream','original_packets','removed_packets', ...
    'retained_packets','windows','mean_phi_formula_sum','mean_phi_time_average'});
disp(Summary);

visibility='off'; if C.show_figures, visibility='on'; end
fig=figure('Name','RapidPD paper equations (10)-(19)','Visible',visibility, ...
    'Color','w','Position',[80 100 1400 520]);
styles={'-','--'};
streamCount=max(arrayfun(@(r)numel(r.metadata.selected_rx_indices),Result));
tiledlayout(1,streamCount,'TileSpacing','compact','Padding','compact');
for s=1:streamCount
    nexttile; hold on;
    for j=1:2
        field=sprintf('stream_%d_phi_formula_sum',s);
        if ismember(field,Result(j).windows.Properties.VariableNames)
            plot(Result(j).windows.stop,Result(j).windows.(field),styles{j}, ...
                'LineWidth',1.4,'DisplayName',Result(j).name);
        end
    end
    grid on; legend('Interpreter','none','Location','best');
    xlabel('Time from recording start (s)'); ylabel('\phi_q = sum_t \psi_n(t)');
    title(sprintf('Rx stream %d: equation (19)',s));
end
sgtitle(sprintf('RapidPD per-Rx scores: %d-layer ACF, no cross-Rx aggregation, lag=%d', ...
    C.paper_acf_layers,C.acf_shift));

out=fullfile(C.path_res,sprintf('rapidpd_paper_%dlayer_%s',C.paper_acf_layers, ...
    char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))));
mkdir(out);
for j=1:2
    writetable(Result(j).filter.table,fullfile(out,sprintf('input_%d_packet_filter.csv',j)));
    writetable(Result(j).windows,fullfile(out,sprintf('input_%d_motion_scores.csv',j)));
end
writetable(Summary,fullfile(out,'summary.csv'));
save(fullfile(out,'rapidpd_paper_results.mat'),'Result','Summary','C','-v7.3');
exportgraphics(fig,fullfile(out,'rapidpd_paper_motion_scores.png'),'Resolution',170);
fprintf('Output: %s\n',out);
if ~C.show_figures, close(fig); end
end

function value=safeMean(x)
if isempty(x), value=NaN; else, value=mean(x,'omitnan'); end
end
