function testRapidPD
% 数学定义、数据质量、因果判决及真实解析器的回归检查。
codeDir=fileparts(fileparts(mfilename('fullpath')));
addpath(codeDir,fullfile(codeDir,'my_function'));
C=configLoad; rng(20260908);
x=randn(5,37); [r,L]=multiLayerACF(x,3);
manual=x;
for layer=1:3
    for row=1:5
        a=conv(manual(row,:),fliplr(manual(row,:)));
        manual(row,:)=a(37:end)/a(37);
    end
    assert(max(abs(L{layer}(:)-manual(:)))<1e-12);
end
assert(max(abs(r(:)-manual(:)))<1e-12);
[~,literal]=multiLayerACF(x,3,'literal');
assert(isequal(literal{1},literal{3}));
assert(all(multiLayerACF(zeros(3,20),3)==0,'all'));
base=1+rand(1,128); csi=repmat(reshape(base,1,1,1,[]),1,2,20,1);
d=winProcess(csi,C); assert(d.score==0 && isnan(d.present));
a=abs(randn(20,128))+1; scale=exp(randn(20,1));
assert(max(abs(csiComplexNorm(a)-csiComplexNorm(a.*scale)),[],'all')<1e-14);
C.threshold=0.2; C.time_aggregation='sum'; C.stream_aggregation='sum';
d=winProcess(reshape(a,1,1,20,128),C);
assert(abs(d.score-d.score_sum)<1e-12);
C=configLoad;
w=winSplit((0:39)'/20,C); assert(numel(w)==2 && all([w.valid]));
t=(0:199)'/20+0.001*sin((0:199)');
w=winSplit(t,C); assert(all([w.valid]));
w=winSplit((0:8)'/20,C); assert(isempty(w));
t=(0:39)'/20; t(t>=0.3 & t<=0.6)=[];
w=winSplit(t,C); assert(~w(1).valid && w(2).valid);
assert(isequaln(labelFilt([1 0 1 0 NaN 1 1 1],3),[NaN;NaN;1;0;NaN;NaN;NaN;1]));
T=table(["a";"b";"c";"d"],ones(4,1),repmat("profile",4,1), ...
    [0.1;0.9;0.2;0.8],[0;1;0;1],["calibration";"calibration";"test";"test"], ...
    'VariableNames',{'file','group','profile','score','label','split'});
M=calibrateThreshold(T); assert(M.test_accuracy==1 && abs(M.threshold-0.5)<1e-12);
T.file(3)="a"; rejected=false;
try, calibrateThreshold(T); catch ME, rejected=strcmp(ME.identifier,'RapidPD:Leakage'); end
assert(rejected);
% Known raw IQ and official format mapping on first supplied record.
root=fileparts(codeDir); sample=fullfile(root,'data','rx_2_260904_154612.csi');
if isfile(sample)
    setupParser(C); [G,A]=fileLoad(sample);
    assert(A.frames==2285 && A.corrupt==0 && numel(G)==1);
    assert(isequal(size(G.csi_raw),[1 2 2285 980]));
    assert(G.csi_raw(1,1,1,1)==complex(single(43),single(29)));
    assert(abs(G.median_rate_hz-20)<0.05);
    assert(numel(G.subcarrier_index)==980 && all(diff(G.subcarrier_index)>0));
end
fprintf('PASS: RapidPD mathematical, quality, smoothing and real-data checks.\n');
end
