function filtered = labelFilt(raw,m)
%LABELFILT 因果多数投票；前 m-1 个窗口为 NaN；无效窗重置连续投票。
raw=raw(:); filtered=nan(size(raw));
validateattributes(m,{'numeric'},{'scalar','integer','positive'});
assert(mod(m,2)==1,'RapidPD:Vote','m 必须为奇数。');
for k=m:numel(raw)
    recent=raw(k-m+1:k);
    if all(isfinite(recent)), filtered(k)=double(sum(recent)>=ceil(m/2)); end
end
end
