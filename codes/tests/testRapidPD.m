function testRapidPD
% 固定20包窗口、Desay FACF 数学路径和论文幅度补偿的回归检查。
codeDir=fileparts(fileparts(mfilename('fullpath')));
addpath(codeDir,fullfile(codeDir,'my_function'));
C=configLoad; rng(20260909);

w=winSplit((0:39)'/20,C);
assert(numel(w)==2 && isequal(w(1).indices,1:20) && isequal(w(2).indices,21:40));
w=winSplit((0:44)'/20,C); assert(numel(w)==2);
w=winSplit((0:18)'/20,C); assert(isempty(w));

a=abs(randn(20,31))+1; scale=exp(randn(20,1));
assert(max(abs(csiComplexNorm(a)-csiComplexNorm(a.*scale)),[],'all')<1e-14);

H=csiComplexNorm(a); got=desayFACF(H,3,1); expected=referenceFACF(H,3,1);
assert(abs(got-expected)<1e-11);
csi=reshape(a,1,1,20,31); detail=winProcess(csi,C);
assert(abs(detail.score-expected)<1e-11);
assert(detail.agc_compensation=="paper_normalization" && isnan(detail.present));

csi2=repmat(csi,1,2,1,1); dMean=winProcess(csi2,C);
C.stream_aggregation='sum'; dSum=winProcess(csi2,C);
assert(abs(dSum.score-2*dMean.score)<1e-11);

fprintf('PASS: fixed 20-packet windows, Desay FACF and paper amplitude compensation.\n');
end

function psi=referenceFACF(H,dmin,shift)
[T,F]=size(H); HH=[H;H]; first=zeros(T*(T-dmin),F-1); row=0;
for p=1:T
    D=H(p,:)-HH(p+(dmin:T-1),:);
    for j=1:size(D,1)
        row=row+1; first(row,:)=referencePositiveACF(D(j,:));
    end
end
first=first-mean(first(:)); second=zeros(size(first,1),F-2);
for j=1:size(first,1), second(j,:)=referencePositiveACF(first(j,:)); end
DD=second-mean(second(:)); DS=circshift(DD,[0 shift]);
psi=sum(DD(:).*DS(:))/sum(DD(:).^2);
end

function y=referencePositiveACF(x)
K=numel(x); y=zeros(1,K-1); zero=sum(x.^2)/K;
for lag=1:K-1
    y(lag)=(sum(x(1:K-lag).*x(1+lag:K))/(K-lag))/zero;
end
end
