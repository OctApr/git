function [rho,layers] = rapidPDMultiLayerACF(HD,n,mode)
%RAPIDPDMULTILAYERACF Implement the visible content of equations (14)-(17).
% First calculate sample self-covariance from HD using (16), then normalize
% it by the zero-lag value using (14). For n>=2, mode='literal' follows the
% printed (17), while mode='recursive' follows the prose/Fig. 5 by applying
% another sample self-covariance operation to the preceding rho sequence.
if nargin<3, mode='recursive'; end
validateattributes(HD,{'numeric'},{'2d','real','finite','nonempty'});
validateattributes(n,{'numeric'},{'scalar','integer','positive'});
assert(size(HD,2)>=2,'RapidPD:Tones','At least two subcarriers are required.');
assert(ismember(string(mode),["recursive","literal"]),'RapidPD:ACFMode','Unknown ACF mode.');

layers=cell(1,n);

% Equations (16) and (14): HD -> sample self-covariance -> rho_1.
gamma=sampleSelfCovariance(double(HD));
rho=normalizeZeroLag(gamma);
layers{1}=rho;

for layer=2:n
    if strcmp(string(mode),"recursive")
        % Candidate interpretation required for distinct curves in Fig. 5.
        gamma=sampleSelfCovariance(rho);
        rho=normalizeZeroLag(gamma);
    else
        % Printed (17): rho_n(t,v)=rho_(n-1)(t,v)/rho_(n-1)(t,0).
        rho=normalizeZeroLag(rho);
    end
    layers{layer}=rho;
end
end

function gamma=sampleSelfCovariance(x)
% Equation (16), apart from the common 1/K factor. Zero-padded FFT gives
% exactly sum_{i=1+k}^K x(i-k)*x(i) for every non-negative linear lag k.
K=size(x,2); N=2^nextpow2(2*K-1);
spectrum=fft(x,N,2);
gamma=real(ifft(spectrum.*conj(spectrum),[],2));
gamma=gamma(:,1:K);
% Equation (16) divides every lag by the same K. It is omitted because it
% cancels exactly in equation (14), gamma(t,v)/gamma(t,0).
end

function rho=normalizeZeroLag(gamma)
denominator=gamma(:,1);
rho=zeros(size(gamma));
valid=denominator>0;
rho(valid,:)=gamma(valid,:)./denominator(valid);
rho(valid,1)=1;
end
