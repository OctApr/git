function T = buildLabeledWindows(resultsFile,labelsFile,outputFile)
% 将人工时间段标签匹配到完整包含于区间的有效窗口。
S=load(resultsFile,'Results');
L=readtable(labelsFile,'TextType','string','Encoding','UTF-8','Delimiter',',');
required={'file','start_s','end_s','label','split'};
assert(all(ismember(required,L.Properties.VariableNames)),'RapidPD:Labels','时间段标签列缺失。');
L.file=replace(string(L.file),'\','/');
for key=["start_s","end_s","label"]
    if ~isnumeric(L.(key)), L.(key)=str2double(string(L.(key))); end
end
L=L(isfinite(L.label),:);
assert(all(ismember(L.label,[0 1])) && all(isfinite(L.start_s)) && ...
    all(isfinite(L.end_s)) && all(L.end_s>L.start_s), ...
    'RapidPD:Labels','已标注区间需要有效起止时间和 0/1 标签。');
assert(all(ismember(string(L.split),["calibration","test"])), ...
    'RapidPD:Labels','已标注区间需要 calibration/test 划分。');
for name=unique(L.file).'
    rows=L(L.file==name,:); rows=sortrows(rows,'start_s');
    assert(all(rows.start_s(2:end)>=rows.end_s(1:end-1)), ...
        'RapidPD:Labels','同一文件的标签时间段不得重叠。');
end
T=table('Size',[0 8],'VariableTypes',{'string','double','string','double','double','string','double','double'}, ...
    'VariableNames',{'file','group','profile','score','label','split','start','stop'});
for fi=1:numel(S.Results)
    R=S.Results{fi}; name=replace(string(R.file),'\','/'); rows=L(L.file==name,:);
    for gi=1:numel(R.groups)
        W=R.groups{gi}.windows;
        for li=1:height(rows)
            keep=W.valid & isfinite(W.score) & W.start>=rows.start_s(li) & W.stop<=rows.end_s(li);
            n=sum(keep);
            if n==0, continue; end
            T=[T;table(repmat(name,n,1),repmat(gi,n,1),W.profile(keep),W.score(keep), ...
                repmat(rows.label(li),n,1),repmat(string(rows.split(li)),n,1),W.start(keep),W.stop(keep), ...
                'VariableNames',T.Properties.VariableNames)]; %#ok<AGROW>
        end
    end
end
if nargin>2, writetable(T,outputFile); end
if isempty(T), warning('RapidPD:NoLabels','没有完整匹配的已标注窗口，尚不能校准。'); end
end
