function testRapidPD
% 固定20包窗口、Desay FACF 数学路径和论文幅度补偿的回归检查。
codeDir=fileparts(fileparts(mfilename('fullpath')));
addpath(codeDir,fullfile(codeDir,'my_function'));
C=configLoad; rng(20260909);
assert(C.paper_acf_layers==3 && configLoad('paper_acf_layers',1).paper_acf_layers==1);
assert(string(configLoad('paper_time_aggregation','mean').paper_time_aggregation)=="mean");

w=winSplit((0:39)'/20,C);
assert(numel(w)==2 && isequal(w(1).indices,1:20) && isequal(w(2).indices,21:40));
w=winSplit((0:44)'/20,C); assert(numel(w)==2);
w=winSplit((0:18)'/20,C); assert(isempty(w));
keep=true(60,1); keep(20)=false;
w=winSplitClean((0:59)'/20,keep,C);
assert(numel(w)==2 && isequal(w(1).indices,21:40) && isequal(w(2).indices,41:60));
assert(~any([w.indices]==20));

a=abs(randn(20,31))+1; scale=exp(randn(20,1));
assert(max(abs(csiComplexNorm(a)-csiComplexNorm(a.*scale)),[],'all')<1e-14);

H=csiComplexNorm(a); got=desayFACF(H,3,1); expected=referenceFACF(H,3,1);
assert(abs(got-expected)<1e-11);
got4=desayFACF(H,3,4); expected4=referenceFACF(H,3,4);
assert(abs(got4-expected4)<1e-11);
csi=reshape(a,1,1,20,31); detail=winProcess(csi,C);
assert(abs(detail.score-expected)<1e-11);
assert(detail.agc_compensation=="paper_normalization" && isnan(detail.present));

csi2=repmat(csi,1,2,1,1); dMean=winProcess(csi2,C);
C.stream_aggregation='sum'; dSum=winProcess(csi2,C);
assert(abs(dSum.score-2*dMean.score)<1e-11);

% Common packet gain must not trigger the spectral-shape filter, while a
% minority of packets with a different frequency profile must be detected.
N=100; F=31; base=1+linspace(0,1,F); synthetic=zeros(1,1,N,F);
for n=1:N
    synthetic(1,1,n,:)=exp(.2*randn)*(base.*(1+.005*randn(1,F)));
end
changed=11:15; alternate=fliplr(base);
for n=changed, synthetic(1,1,n,:)=exp(.2*randn)*(alternate.*(1+.005*randn(1,F))); end
filter=packetSpectralOutliers(synthetic);
assert(all(filter.removed(changed)) && sum(filter.removed)<=numel(changed)+1);

fprintf('PASS: fixed 20-packet windows, Desay FACF and paper amplitude compensation.\n');

% Paper equations (10)-(19): the residual must have zero temporal mean for
% every tone, and the reported stream statistic is the packet mean at lag.
paperInput=reshape(abs(randn(1,2,20,31))+0.1,1,2,20,31);
paper=rapidPDPaperWindow(paperInput,'acf_layers',3,'acf_lag',1);
for s=1:numel(paper.stream_details)
    residual=paper.stream_details{s}.residual;
    assert(max(abs(mean(residual,1)))<1e-14);
    assert(abs(paper.stream_phi(s)-mean(paper.packet_statistics(:,s)))<1e-14);
    assert(abs(paper.stream_phi_formula_sum(s)-sum(paper.packet_statistics(:,s)))<1e-14);
end
assert(abs(paper.score-sum(paper.stream_phi_formula_sum))<1e-14);
assert(abs(paper.Phi-paper.Phi_formula_sum)<1e-14);
assert(abs(paper.Phi_formula_sum-20*paper.Phi_text_average)<1e-12);

% FFT implementation must equal the explicit product sum in equation (16)
% after the common 1/K factor cancels through equation (14).
smallHD=randn(4,9);
[rhoFFT,~]=rapidPDMultiLayerACF(smallHD,1,'literal');
rhoDirect=zeros(size(smallHD));
for t=1:size(smallHD,1)
    for k=0:size(smallHD,2)-1
        gamma=sum(smallHD(t,1:end-k).*smallHD(t,1+k:end))/size(smallHD,2);
        gamma0=sum(smallHD(t,:).^2)/size(smallHD,2);
        rhoDirect(t,k+1)=gamma/gamma0;
    end
end
assert(max(abs(rhoFFT-rhoDirect),[],'all')<1e-12);
literal=rapidPDPaperWindow(paperInput,'acf_layers',3,'acf_lag',1,'acf_mode','literal');
oneLayer=rapidPDPaperWindow(paperInput,'acf_layers',1,'acf_lag',1,'acf_mode','literal');
assert(max(abs(literal.packet_statistics-oneLayer.packet_statistics),[],'all')<1e-14);
fprintf('PASS: RapidPD paper normalization, window benchmark and three-layer ACF.\n');
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
