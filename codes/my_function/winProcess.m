function Detail = winProcess(csi_split,Config)
%WINPROCESS Desay：论文式幅度补偿->包间差分->两层单边 ACF->整体相关。
ntx=size(csi_split,1); nrx=size(csi_split,2);
T=size(csi_split,3); F=size(csi_split,4);
assert(T==Config.window_packets,'RapidPD:Window', ...
    '每个窗口必须严格包含 %d 个包，当前为 %d。',Config.window_packets,T);
assert(F>=3,'RapidPD:Window','至少需要三个子载波。');
stream=nan(1,ntx*nrx);
for rx=1:nrx
    for tx=1:ntx
        s=(rx-1)*ntx+tx;
        H=csiComplexNorm(reshape(csi_split(tx,rx,:,:),T,F));
        stream(s)=desayFACF(H,Config.pair_min_distance,Config.acf_shift);
    end
end
Detail=struct('packet_statistics',[],'stream_statistics',stream, ...
    'score_mean',mean(stream),'score_sum',sum(stream), ...
    'score_stream_sum',sum(stream),'score',NaN,'present',NaN, ...
    'agc_compensation','paper_normalization');
if strcmp(Config.stream_aggregation,'sum'), Detail.score=sum(stream);
else, Detail.score=mean(stream); end
if isfinite(Config.threshold), Detail.present=double(Detail.score>=Config.threshold); end
end
