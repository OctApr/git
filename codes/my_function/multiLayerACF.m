function [rho, layers] = multiLayerACF(D,n,mode)
%MULTILAYERACF 子载波维的线性、biased 自相关；保留零延迟项。
% 式(16)固定分母 K 与按零延迟归一化约掉。FFT 零填充避免循环相关。
% recursive: 依据文字每层重新自相关；literal: 直接照式(17)只归一化。
if nargin<3, mode='recursive'; end
validateattributes(D,{'numeric'},{'2d','real','finite','nonempty'});
validateattributes(n,{'numeric'},{'scalar','integer','positive'});
assert(size(D,2)>=2,'RapidPD:Tones','至少需要两个子载波。');
assert(ismember(string(mode),["recursive","literal"]),'RapidPD:ACFMode','模式无效。');
rho=double(D); layers=cell(1,n);
for k=1:n
    if k==1 || strcmp(mode,'recursive')
        K=size(rho,2); N=2^nextpow2(2*K-1);
        spectrum=fft(rho,N,2);
        a=real(ifft(spectrum.*conj(spectrum),[],2)); a=a(:,1:K);
        energy=a(:,1);
        rho=zeros(size(a)); good=energy>0;
        rho(good,:)=a(good,:)./energy(good);
        rho(good,1)=1;
    end
    layers{k}=rho;
end
end
