function rx = selectRxIndices(nrx,Config)
%SELECTRXINDICES 返回本次参与算法的接收链下标。
if isempty(Config.rx_indices)
    rx=1:nrx;
else
    rx=double(Config.rx_indices(:).');
    assert(all(rx<=nrx),'RapidPD:RxIndex', ...
        '数据只有 %d 条接收链，但 rx_indices 包含 %d。',nrx,max(rx));
end
end
