function Model = calibrateThreshold(scoreTable,outputFile)
%CALIBRATETHRESHOLD 仅对标注为 calibration 的独立数据选择阈值。
% 输入 CSV/table 必须有 file,group,profile,score,label,split 六列。
% 同一文件不能同时用于 calibration/test；相同物理配置须分别校准。
% label=0/1，split="calibration"/"test"；未知/切换窗口不要加入。
if ~istable(scoreTable)
    scoreTable=readtable(scoreTable,'TextType','string','Encoding','UTF-8','Delimiter',',');
end
assert(all(ismember({'file','group','profile','score','label','split'},scoreTable.Properties.VariableNames)), ...
    'RapidPD:Labels','需要 file,group,profile,score,label,split 列。');
T=scoreTable;
assert(numel(unique(string(T.profile)))==1,'RapidPD:Profile', ...
    '不同发射端或物理配置应分别校准，不可混用阈值。');
assert(all(ismember(T.label,[0 1])) && all(isfinite(T.score)),'RapidPD:Labels','分数和标签无效。');
cal=string(T.split)=="calibration"; test=string(T.split)=="test";
assert(all(cal|test),'RapidPD:Labels','split 只能为 calibration/test。');
assert(isempty(intersect(string(T.file(cal)),string(T.file(test)))), ...
    'RapidPD:Leakage','同一文件不可跨校准集和测试集。');
assert(numel(unique(T.label(cal)))==2,'RapidPD:Labels','校准集需要有人和无人两类。');
s=T.score(cal); y=T.label(cal); u=unique(s);
candidates=[u(1)-max(eps(u(1)),1e-12);(u(1:end-1)+u(2:end))/2;u(end)+max(eps(u(end)),1e-12)];
ba=zeros(size(candidates));
for k=1:numel(candidates)
    p=s>=candidates(k); ba(k)=(mean(p(y==1))+mean(~p(y==0)))/2;
end
[best,k]=max(ba);
Model=struct('threshold',candidates(k),'calibration_balanced_accuracy',best, ...
    'profile',string(T.profile(1)), ...
    'calibration_windows',sum(cal),'test_windows',sum(test), ...
    'test_accuracy',NaN,'test_tpr',NaN,'test_tnr',NaN);
if any(test)
    p=T.score(test)>=Model.threshold; y=T.label(test);
    Model.test_accuracy=mean(p==y); Model.test_tpr=mean(p(y==1)); Model.test_tnr=mean(~p(y==0));
end
if nargin>1, save(outputFile,'Model'); end
end
