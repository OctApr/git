function Detail = rapidPDPaperWindow(csi_split,varargin)
%RAPIDPDPAPERWINDOW RapidPD equations (10)-(19) for one 20-packet window.
% Order: amplitude normalization -> per-tone window mean -> residual ->
% three recursive one-sided ACF layers -> non-zero-lag packet statistic.
p=inputParser;
addParameter(p,'acf_layers',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
addParameter(p,'acf_lag',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
parse(p,varargin{:});

ntx=size(csi_split,1); nrx=size(csi_split,2);
T=size(csi_split,3); F=size(csi_split,4);
assert(T==20,'RapidPD:Window','Paper experiment requires exactly 20 packets.');
assert(p.Results.acf_lag<F,'RapidPD:Lag','ACF lag must be smaller than tone count.');

S=ntx*nrx; psi=nan(T,S); phi=nan(1,S);
details=cell(1,S);
for rx=1:nrx
    for tx=1:ntx
        s=(rx-1)*ntx+tx;
        H=csiComplexNorm(reshape(csi_split(tx,rx,:,:),T,F)); % (10)-(11)
        Hbar=mean(H,1);                                     % (12)
        HD=H-Hbar;                                          % (13)
        if max(abs(HD(:)))<=32*eps(max(H(:))), HD(:)=0; end
        [rho,layers]=rapidPDMultiLayerACF(HD,p.Results.acf_layers);
        psi(:,s)=rho(:,p.Results.acf_lag+1);                % (18)
        phi(s)=mean(psi(:,s)); % paper calls phi an average; (19) omits 1/T
        details{s}=struct('normalized_amplitude',H,'window_mean',Hbar, ...
            'residual',HD,'acf_layers',{layers});
    end
end

% The paper sums Tx-Rx streams for overall Phi.  Keep a mean-scale score as
% well because AX210 has 1Tx-2Rx rather than the paper's 2Tx-1Rx topology.
Detail=struct('packet_statistics',psi,'stream_phi',phi, ...
    'score',sum(phi),'score_stream_sum',sum(phi), ...
    'score_stream_mean',mean(phi),'stream_details',{details}, ...
    'acf_layers',p.Results.acf_layers,'acf_lag',p.Results.acf_lag, ...
    'time_aggregation','mean','stream_aggregation','sum');
end
