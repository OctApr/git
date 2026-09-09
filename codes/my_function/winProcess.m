function Detail = winProcess(csi_split,Config)
%WINPROCESS 论文 III-B/C/D：归一化->窗口背景->多层 ACF->lag1。
ntx=size(csi_split,1); nrx=size(csi_split,2);
T=size(csi_split,3); F=size(csi_split,4);
assert(T>=2 && F>=2,'RapidPD:Window','窗口至少需要两个包和两个子载波。');
psi=zeros(T,ntx*nrx); stream=nan(1,ntx*nrx);
for rx=1:nrx
    for tx=1:ntx
        s=(rx-1)*ntx+tx;
        H=csiComplexNorm(reshape(csi_split(tx,rx,:,:),T,F));
        D=H-mean(H,1); % 式(12)-(13)，仅当前窗口，不使用未来窗口。
        % Identical normalized rows can leave roundoff after temporal mean.
        if max(abs(D(:)))<=32*eps(max(H(:))), D(:)=0; end
        rho=multiLayerACF(D,Config.acf_layers,Config.acf_mode);
        psi(:,s)=rho(:,2); % 式(18) 一个子载波索引延迟，不取绝对值。
        if strcmp(Config.time_aggregation,'sum'), stream(s)=sum(psi(:,s));
        else, stream(s)=mean(psi(:,s)); end
    end
end
Detail=struct('packet_statistics',psi,'stream_statistics',stream, ...
    'score_mean',mean(psi(:)),'score_sum',sum(psi(:)), ...
    'score_stream_sum',sum(mean(psi,1)),'score',NaN,'present',NaN);
if strcmp(Config.stream_aggregation,'sum'), Detail.score=sum(stream);
else, Detail.score=mean(stream); end
if isfinite(Config.threshold), Detail.present=double(Detail.score>=Config.threshold); end
end
