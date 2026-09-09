# 第二轮 Desay FACF 复现约定

参考 RapidPD 论文和 DesayCPD 提交 `f705e2feb9ece2d54b0bc98d8dd68c6dbfea644a`。本轮以 Desay 的有效代码路径为算法主体，只在 AX210 无法提供等价 AGC 字段时采用论文幅度补偿。

| 环节 | 第二轮实现 |
|---|---|
| CSI 组织 | `[Tx,Rx,Packet,Subcarrier]` |
| 分窗 | 连续20包一窗，步长20包，末尾不足20包丢弃 |
| 幅度补偿 | 每条流、每个包的幅值除以该包所有有效子载波幅值和 |
| 环境变化 | 当前包与循环窗口中距离至少3包的其他包两两差分 |
| 第一层 ACF | `xcorr(...,'unbiased')` 等价计算，只保留正延迟，除以零延迟 |
| 层间处理 | 第一层矩阵减去全矩阵标量均值 |
| 第二层 ACF | 再次计算单边正延迟 unbiased ACF，除以零延迟 |
| 时间汇总 | 第二层矩阵全局去均值后循环右移1列，计算全部元素的归一化内积 |
| 跨流汇总 | 默认平均，与 Desay `judgFun` 一致；可显式配置为求和 |
| 判决 | 默认阈值为 NaN，只输出曲线；阈值必须用当前数据重新校准 |

设单流归一化幅度矩阵为 `H`，维度为 `20 × F`。Desay 使用 `HH=[H;H]`，对每个包 `i` 构造：

```matlab
D = H(i,:) - HH(i+(3:19),:);
```

因此每条流每窗形成 `20 × 17 = 340` 个包对差分。窗口尾部会与窗口开头形成包对，这是 Desay 源码的循环边界。

长度为 `K` 的序列完整互相关长度为 `2K-1`。代码只保留 `lag=1...K-1`，不保留负延迟和零延迟，所以第一层长度为 `F-1`，第二层长度为 `F-2`。两层使用线性 ACF 和 `unbiased` 分母；最终 `circshift` 使用循环边界。

最终单流分数为：

```matlab
DD = second - mean(second(:));
DS = circshift(DD,[0 1]);
psi = sum(DD(:).*DS(:)) / sum(DD(:).^2);
```

这一步已经把窗口内所有包对和频率延迟点整体压缩成一个值，因此没有额外的逐时间平均或求和。

## AX210 的 AGC 边界

Desay 数据含每条接收链的 `agc`/`sts_tot_gain`，并执行：

```matlab
csi(:,i,:) = csi(:,i,:) / 10^(agc(i)/20);
```

当前 PicoScenes AX210 帧的 `MVMExtra` 只有 `IQDataSize`、`FTMClock`、`NumTones`、`RateNFlags`；`RxSBasic` 提供 RSSI 和噪声，但没有等价的 AGC 总增益字段。RSSI 不能在没有硬件标定关系的情况下直接替代 AGC。

所以本实现采用 RapidPD 论文式幅度补偿：每包、每流除以所有有效子载波幅值和。它能严格抵消同一包同一流的共同乘性增益，但不声称已经复现 Desay 的逐接收链 AGC 逆补偿。输出元数据中记录：

```text
native_agc_available = false
amplitude_compensation = paper_normalization
```

## 与第一轮的隔离

第一轮代码保存在 Git 标签 `rapidpd-round1`。第二轮位于分支 `codex/rapidpd-round2`。两轮分数定义和阈值尺度不同，不能混用第一轮结果或阈值。
