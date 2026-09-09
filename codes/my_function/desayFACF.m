function [psi, Layers] = desayFACF(H,dmin,shift)
%DESAYFACF 复现 Desay winProcessFACF 的有效计算路径。
% H: [20 packets, F tones]；输出为一个流在该窗口内的整体归一化相关。
validateattributes(H,{'numeric'},{'2d','real','finite','nonempty'});
[T,F]=size(H);
assert(T>dmin && F>=3,'RapidPD:FACF','数据包或子载波数量不足。');
assert(shift<size(H,2)-2,'RapidPD:FACF','循环移位必须小于第二层 ACF 长度。');

HH=[H;H];
first=zeros(T*(T-dmin),F-1);
row=0;
for iPkt=1:T
    D=H(iPkt,:)-HH(iPkt+(dmin:T-1),:);
    rows=row+(1:size(D,1));
    first(rows,:)=positiveUnbiasedACF(D);
    row=row+size(D,1);
end
first=first-mean(first(:));

second=positiveUnbiasedACF(first);
DD=second-mean(second(:));
DS=circshift(DD,[0 shift]);
den=sum(DD(:).^2);
if den==0
    psi=NaN;
else
    psi=sum(DD(:).*DS(:))/den;
end
Layers=struct('first',first,'second',second,'centered_second',DD);
end

function out = positiveUnbiasedACF(x)
% 等价于 xcorr(x,'unbiased') 的正延迟部分除以零延迟。
x=double(x); K=size(x,2); N=2^nextpow2(2*K-1);
spectrum=fft(x,N,2);
a=real(ifft(spectrum.*conj(spectrum),[],2));
unbiased=a(:,1:K)./(K-(0:K-1));
zero=unbiased(:,1);
out=unbiased(:,2:end)./zero;
out(zero==0,:)=NaN;
end
