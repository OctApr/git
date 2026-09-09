function H = csiComplexNorm(csi)
%CSICOMPLEXNORM 论文式(10)-(11)：各包幅值除以子载波幅值和。
A=abs(double(csi));
assert(ismatrix(A) && all(isfinite(A(:))),'RapidPD:InvalidCSI','CSI 含非有限值。');
s=sum(A,2);
assert(all(s>0),'RapidPD:InvalidCSI','全零 CSI 包不能归一化。');
H=A./s;
end
