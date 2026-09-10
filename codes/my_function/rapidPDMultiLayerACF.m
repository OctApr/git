function [rho,layers] = rapidPDMultiLayerACF(HD,n)
%RAPIDPDMULTILAYERACF Paper-visible recursive subcarrier-domain ACF.
% Each row is one packet.  Keep the non-negative linear lags (including
% lag zero), normalize by lag zero, and feed the result to the next layer.
validateattributes(HD,{'numeric'},{'2d','real','finite','nonempty'});
validateattributes(n,{'numeric'},{'scalar','integer','positive'});
assert(size(HD,2)>=2,'RapidPD:Tones','At least two subcarriers are required.');

rho=double(HD); layers=cell(1,n);
for layer=1:n
    K=size(rho,2); N=2^nextpow2(2*K-1);
    spectrum=fft(rho,N,2);
    acf=real(ifft(spectrum.*conj(spectrum),[],2));
    acf=acf(:,1:K);                 % lag 0,...,K-1 (one-sided)
    energy=acf(:,1);
    next=zeros(size(acf));
    valid=energy>0;
    next(valid,:)=acf(valid,:)./energy(valid);
    next(valid,1)=1;
    rho=next;
    layers{layer}=rho;
end
end
