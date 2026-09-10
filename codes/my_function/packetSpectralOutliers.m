function Filter = packetSpectralOutliers(csi_raw,varargin)
%PACKETSPECTRALOUTLIERS Detect packet-wide normalized spectral-shape changes.
% The score is insensitive to a common multiplicative gain because every
% Tx/Rx stream is first normalized by its sum over subcarriers.
p=inputParser;
addParameter(p,'mad_multiplier',6,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(p,varargin{:});

ntx=size(csi_raw,1); nrx=size(csi_raw,2);
N=size(csi_raw,3); F=size(csi_raw,4); S=ntx*nrx;
assert(N>=3 && F>=3,'RapidPD:OutlierFilter','CSI dimensions are too small.');
streamError=nan(N,S); s=0;
for rx=1:nrx
    for tx=1:ntx
        s=s+1;
        H=abs(reshape(csi_raw(tx,rx,:,:),N,F));
        rowSum=sum(H,2);
        assert(all(isfinite(rowSum)&rowSum>0),'RapidPD:OutlierFilter', ...
            'CSI contains an invalid all-zero or non-finite packet.');
        H=H./rowSum;
        reference=median(H,1);
        reference=reference/sum(reference);
        streamError(:,s)=sqrt(sum((H-reference).^2,2))/norm(reference);
    end
end
score=max(streamError,[],2);
center=median(score);
robustSigma=1.4826*median(abs(score-center));
threshold=center+p.Results.mad_multiplier*robustSigma;
removed=score>threshold;
Filter=struct('keep',~removed,'removed',removed,'score',score, ...
    'stream_error',streamError,'center',center,'robust_sigma',robustSigma, ...
    'threshold',threshold,'mad_multiplier',p.Results.mad_multiplier, ...
    'table',table((1:N).',score,removed,'VariableNames', ...
    {'packet_index','spectral_shape_error','removed'}));
end
