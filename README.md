# AX210 CSI 存在检测：RapidPD MATLAB 实现

读取两台 AX210 主机一发一收采集的 PicoScenes `.csi`。第二轮算法复现 Desay FACF 的固定20包、循环包对差分、两层单边 ACF 和整体归一化相关；AX210 没有输出 Desay 所需的 AGC 增益字段，因此幅度补偿采用 RapidPD 论文的逐包幅值和归一化。

## 快速运行

### 推荐：每次只处理两份工作区数据

先将两个 `.csi` 按现有 PicoScenes 导入方式拖拽解析到 MATLAB 工作区，再传入**实际变量名**：

```matlab
addpath('C:/Users/13059/Documents/ChatGPT/git/codes');
R = runPair(rx_2_260904_171616, rx_2_260904_180803);
```

当前数据实际为1个发送空间流、2个接收链。先只使用Rx1时：

```matlab
R = runPair(rx_2_260904_171616, rx_2_260904_180803, 'rx_indices', 1);
```

`rx_indices=[]`（默认）分别计算两条接收链后取平均；`rx_indices=2`只使用Rx2。数据不会沿天线维拼接成一条子载波序列。

对PicoScenes默认插值后的1001-tone工作区变量，可运行频率移位量对照：

```matlab
[Summary,Result,fig,outputDir] = compareShifts( ...
    rx_2_260904_171616_1M, rx_2_260904_180803_NO);
```

该实验同时计算`shift=1`（78.125 kHz）和`shift=4`（312.5 kHz）。它要求输入为80 MHz HE、`SubcarrierIndex=-500:500`的数据，避免把有缺口的980-tone序列位置误当成均匀物理频率间隔。

只处理这两个输入，不扫描 `data`，不重新读取工作区变量对应的源文件。图会立即显示，`R(1).windows` 和 `R(2).windows` 是两份逐窗口分数表；PNG、CSV、MAT 同时保存到 `results/pair_时间/`。每窗严格使用20个连续包，末尾不足20包的数据丢弃；横轴来自实际时间戳，不做有人/无人判决。

也可以直接传入两个文件路径，省去先解析到工作区的步骤：

```matlab
R = runPair('F:/BaiduNetdiskDownload/csi/csi (1)一人/rx_2_260904_171616.csi', ...
            'F:/BaiduNetdiskDownload/csi/csi(2)二人+无人/rx_2_260904_180803.csi');
```

支持 `read_rxs_log` 的原始 cell、原始 struct 数组，以及 `parseCSIFile/parseRXSBundle` 产生的合并 struct 或 cell 包装。裸 CSI 矩阵缺少时间戳/维度元数据，不接受。若工具箱返回自定义对象，请使用其原始 cell/struct 导入形式。

**工作区导入与路径导入的区别：**官方默认拖拽解析通常包含插值预处理；工作区输入保留这些已处理数值，并在命令行提示。路径输入继续使用本项目的 `RXSParser(...,false)` 关闭插值，因此两种导入方式的数值可能不同；对比两份实验时应使用一致的导入方法。

### 批量处理整个 data 目录

已在本机 MATLAB R2024b 和现有 PicoScenes 官方 MEX 上验证。算法不依赖 Signal Processing Toolbox 或 Parallel Computing Toolbox；`xcorr(...,'unbiased')` 的等价计算由 MATLAB 内置 FFT 实现。

在项目根目录运行：

```matlab
addpath('codes');
[summary, results] = main;
```

默认处理 `data` 下全部 `.csi` 文件，保存到 `results/运行时间/`：每个物理配置分组的窗口 CSV、诊断 PNG、`summary.csv`、`results.mat` 和 `config.json`。默认阈值未设置，只有分数，二值判决为 NaN。

运行完成会显示分数对比图，不需要先设置阈值。默认只选 HE-SU 帧（`packet_formats=3`），避免把周围设备的普通 Wi-Fi 流量混进实验；探索所有帧可传入 `'packet_formats',[]`。`'source_mac','xx:xx:xx:xx:xx:xx'` 可以进一步限定实验发射端。

只想看已有结果，不重跑采集文件：

```matlab
addpath('codes');
plotComparison;                 % 自动打开最新批处理的分数对比图
plotProcessing('data/csi(2)二人+无人/rx_2_260904_180803.csi');
```

无标签对比图使用文件名作为图例，不自动赋予有人/无人含义。`plotProcessing` 显示原始幅值、论文式幅度补偿、包对差分、第一层单边 ACF、第二层单边 ACF和逐窗口整体相关分数，导出到 `results/processing/`。全部是真实数据，不为贴近论文图而调整数值。

当前按用户指定，结果中同时存在 `171616` 和 `180803` 时，`plotComparison` 默认比较这两份文件；也可直接运行 `compareRecordings`。两条线各保留原有窗口数，不截短或补点。两份数据的源 MAC 不同，图只用于观察原始处理分数，不代表控制了全部采集条件的分类验证。其他数据集仍使用上述自动选择规则。

批处理不希望弹出图窗时用 `main('show_figures',false)`；连图文件也不需要时用 `main('make_plots',false)`。

```matlab
% 先处理单个样例
[summary, results] = main('file_pattern','rx_2_260904_171616.csi');

% 显式传入探索性阈值；0.43 来自论文，未在当前数据上校准
main('threshold',0.43,'file_pattern','rx_2_260904_171616.csi');

% 默认跨流平均与 Desay 一致；也可显式查看跨流求和尺度
main('stream_aggregation','sum');

% 运行回归检查
addpath('codes/tests'); testRapidPD;
```

工具箱未加入路径时：

```matlab
main('parser_path','C:/Users/13059/PicoScenes-MATLAB-Toolbox-Core');
```

其他机器需按 [PicoScenes 官方说明](https://github.com/wifisensing/PicoScenes-MATLAB-Toolbox-Core) 安装/编译 `RXSParser`。本工程不自动修改系统工具箱文件。

## 文件组织

```text
codes/
  main.m                    批处理入口、CSV/MAT/图表输出
  calibrateThreshold.m      有真实标签后的阈值校准
  buildLabeledWindows.m     将人工时间段标签转换为窗口标签
  my_function/
    configLoad.m            所有可配置参数
    setupParser.m           检查官方解析器
    fileLoad.m              逐帧读取、数据审计、物理配置分组
    winSplit.m              连续数据包固定20包分窗
    csiComplexNorm.m        论文式逐包幅值和归一化
    desayFACF.m             包对差分、两层单边ACF与整体相关
    winProcess.m            多天线流单窗分数计算
    labelFilt.m             因果多数投票
    resDisp.m               诊断图
  tests/testRapidPD.m        数学与真实数据检查
data/                       29 个采集文件（仅本地，不提交 Git）
  manifest.csv              路径、字节数、SHA-256
  labels.csv                待补全的时间段标签
  docs/algorithm.md           公式映射、歧义、适配边界
  docs/reproduction_matrix.md 每轮论文/Desay/当前代码复现状态表
results/                    自动生成，不提交 Git
```

## 校准与验证

先补充 `data/labels.csv` 中无人/有人时间段。`buildLabeledWindows` 会选择完全落在已知区间内的有效窗口，转换为供 `calibrateThreshold` 使用的窗口级标注表。

选择同一 profile（发射端/物理参数）、同一算法配置的分数。不同 profile 分别校准；按文件划分 calibration/test，避免重叠窗口或同一采集文件泄漏。函数拒绝混合 profile 和同文件跨数据集。

```matlab
T = buildLabeledWindows('results/运行时间/results.mat','data/labels.csv','labeled_windows.csv');
assert(~isempty(T),'请先补全真实标签');
% 若有多个 profile，先选一个，分别校准
T = T(T.profile == T.profile(1),:);
M = calibrateThreshold(T,'threshold_model.mat');
main('threshold',M.threshold,'threshold_profile',M.profile);
```

该调用仅对匹配 profile 应用阈值，其他组继续只输出分数。算法参数必须与校准时一致。校准目标为原始单窗口的 balanced accuracy；`test_*` 也是未平滑窗口指标，不冒充 3 窗投票指标。只有一类标签时拒绝校准。

`0/1` 是存在与否，不是人数；`NaN` 表示未设置阈值、数据无效或投票未预热，不能解释成无人。

## 与论文的对应

默认流程为论文式逐包幅值和归一化 → Desay 循环包对差分 → 两层正延迟 `unbiased` ACF → 循环移位 → 窗口整体归一化相关 → 跨流平均。完整定义和 AGC 边界见 [算法说明](docs/algorithm.md)。当前未设置分类阈值，曲线不能直接解释为检测准确率。
