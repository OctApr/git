function frames = workspaceFrames(input)
%WORKSPACEFRAMES 支持 read_rxs_log 原始 cell、struct 数组和 parseCSIFile 合并 bundle。
% 只读取输入变量，不读取文件、不扫描 base workspace，不使用 eval。
if iscell(input)
    frames={};
    for i=1:numel(input)
        part=workspaceFrames(input{i}); frames=[frames;part(:)]; %#ok<AGROW>
    end
    return;
end
assert(isstruct(input) && ~isempty(input),'RapidPD:WorkspaceFormat', ...
    '请传入 PicoScenes 的原始 cell/struct 或 parseCSIFile 合并 struct；不接受不带元数据的裸矩阵。');
if numel(input)>1
    frames=workspaceFrames(num2cell(input)); return;
end
assert(all(isfield(input,{'CSI','RxSBasic','StandardHeader'})), ...
    'RapidPD:WorkspaceFormat','变量必须包含 CSI、RxSBasic 和 StandardHeader。');
assert(isfield(input.RxSBasic,'SystemTime'),'RapidPD:WorkspaceFormat','需要 SystemTime 时间戳。');
N=numel(input.RxSBasic.SystemTime);
if N==1
    frames={input}; return;
end
assert(N>1 && isscalar(input.CSI),'RapidPD:WorkspaceFormat','合并数据格式不支持。');
assert(size(input.CSI.CSI,1)==N,'RapidPD:WorkspaceFormat', ...
    '合并 CSI 应为 [包数,展平的CSI]，与官方 parseRXSBundle 格式一致。');
frames=cell(N,1);
for k=1:N
    r=struct;
    for field=["CSI","RxSBasic","StandardHeader"]
        original=input.(field); unpacked=struct;
        keys=fieldnames(original);
        for i=1:numel(keys)
            key=keys{i}; value=original.(key);
            if isnumeric(value) || islogical(value)
                if isempty(value), unpacked.(key)=value;
                elseif size(value,1)==N, unpacked.(key)=value(k,:);
                else
                    error('RapidPD:WorkspaceFormat','字段 %s.%s 的包数不一致。',field,key);
                end
            else
                % Nested unused header metadata is kept without affecting CSI.
                unpacked.(key)=value;
            end
        end
        r.(field)=unpacked;
    end
    frames{k}=r;
end
end
