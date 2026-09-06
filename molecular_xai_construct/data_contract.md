# 数据合同

## 原始基准

- 论文：Jiménez-Luna et al., 2022, DOI `10.1021/acs.jcim.1c01163`。
- 代码：`https://github.com/josejimenezluna/xaibench_tf`。
- commit：`c8a32c4d9253a1531f5ce1fefba47b0509cd21e0`。
- 许可：仓库 `LICENSE` 为 AGPL-3.0。
- 原数据记录：ETH Research Collection handle `20.500.11850/504716`，item UUID `f2e4a104-7590-460c-9c42-fad977f57141`。
- 数据记录：大文件链接位于 `libdrive.ethz.ch`。

## 扩展数据

- 论文：Amara et al., 2023, DOI `10.1186/s13321-023-00733-9`。
- 代码：`https://github.com/microsoft/molucn`。
- commit：`5c9c689c1c9ca1dc12b2c933648bf6e30edf70c3`。
- 数据：Figshare `10.6084/m9.figshare.21215477.v3`，CC BY 4.0。

| 文件 | Figshare file id | 字节数 | MD5 | 用途 |
|---|---:|---:|---|---|
| `data.tar.gz` | 37624043 | 868,218,523 | `cab828a6c6570044e7d9fb9914d38098` | 350 targets 的 pair/ground-truth/MCS 对象 |
| `results.tar.gz` | 37691352 | 12,146,364 | `b83ec75c28c5eb87a1b19ca20bb12507` | 发布的逐 target/method/threshold 分数 |
| `logs.tar.gz` | 37691397 | 265,689 | `ba6e100103955771d1441d983a963d32` | 训练日志 |
| `colors.tar.gz` | 41040914 | 982,773,003 | `1437041174049a9b7f14e47f758847e9` | 预计算 attribution 颜色 |

## 分析范围

- 主分析使用 test split。
- 纳入 MCS 50%–95%，固定 5% 间隔。
- 纳入发布结果中有效的 350 targets；每个终点逐项报告缺失，不做结果驱动填补。
- 比较 GNN explainers 时在同一 loss 内进行；RF 不混入 GNN 内部排名。
- 2023 年扩展包含 350 个 targets，2022 年原始数据包含 738 个 series。
- 数值来自本地计算结果。
