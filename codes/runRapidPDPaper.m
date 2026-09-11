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
    stream_phi=nan(N,size(csi,1)*size(csi,2));
    for wi=1:N
        idx=kept((wi-1)*20+(1:20));
        d=rapidPDPaperWindow(csi(:,:,idx,:), ...
            'acf_layers',3,'acf_lag',C.acf_shift,'acf_mode','recursive');
        literal=rapidPDPaperWindow(csi(:,:,idx,:), ...
            'acf_layers',3,'acf_lag',C.acf_shift,'acf_mode','literal');
        start(wi)=g.time_seconds(idx(1)); stop(wi)=g.time_seconds(idx(end));
        max_gap(wi)=max(diff(g.time_seconds(idx)));
        phi_mean(wi)=d.diagnostic_stream_mean;
        Phi_text_average(wi)=d.Phi_text_average;
        Phi_formula_sum(wi)=d.Phi_formula_sum;
        Phi_literal_formula_sum(wi)=literal.Phi_formula_sum;
        stream_phi(wi,:)=d.stream_phi;
    end
    W=table(start,stop,max_gap,repmat(20,N,1),phi_mean,Phi_text_average, ...
        Phi_formula_sum,Phi_literal_formula_sum, ...
        'VariableNames',{'start','stop','max_gap','packets','diagnostic_stream_mean', ...
        'Phi_text_average','Phi_formula_sum','Phi_literal_formula_sum'});
    for s=1:size(stream_phi,2), W.(sprintf('stream_%d_phi',s))=stream_phi(:,s); end
    metadata=rmfield(g,{'csi_raw','packet_csi'});
    metadata.selected_rx_indices=selectedRx; metadata.audit=audit;
    metadata.processing_order={'remove_spectral_outliers','amplitude_normalize', ...
        'per_tone_window_mean','subtract_window_mean','three_layer_acf','lag_statistic'};
    Result(j)=struct('name',names(j),'filter',filter,'windows',W,'metadata',metadata);
    rows(end+1,:)={names(j),size(csi,3),sum(filter.removed),sum(filter.keep),N, ... %#ok<AGROW>
        safeMean(Phi_text_average),safeMean(Phi_formula_sum),safeMean(Phi_literal_formula_sum)};
    fprintf('%s: removed %d/%d, retained %d, windows %d, recursive Phi %.6f, equation-(19) Phi %.6f, literal-(17) Phi %.6f\n', ...
        names(j),sum(filter.removed),numel(filter.removed),sum(filter.keep),N, ...
        safeMean(Phi_text_average),safeMean(Phi_formula_sum),safeMean(Phi_literal_formula_sum));
end
Summary=cell2table(rows,'VariableNames',{'input','original_packets','removed_packets', ...
    'retained_packets','windows','mean_Phi_text_average','mean_Phi_formula_sum', ...
    'mean_Phi_literal_formula_sum'});
disp(Summary);

visibility='off'; if C.show_figures, visibility='on'; end
fig=figure('Name','RapidPD paper equations (10)-(19)','Visible',visibility, ...
    'Color','w','Position',[100 100 1150 520]); hold on;
styles={'-','--'};
for j=1:2
    plot(Result(j).windows.stop,Result(j).windows.Phi_text_average,styles{j}, ...
        'LineWidth',1.4,'DisplayName',Result(j).name);
end
grid on; legend('Interpreter','none','Location','best');
xlabel('Time from recording start (s)'); ylabel('Overall motion statistic \Phi (stream sum)');
title(sprintf('RapidPD intended path: recursive 3-layer ACF, lag=%d, no division by stream count',C.acf_shift));

out=fullfile(C.path_res,['rapidpd_paper_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
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
