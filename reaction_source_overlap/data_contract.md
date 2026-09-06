# 数据与来源合同

状态：`stopped_source_gate_failed`
冻结时间：2026-08-09（Asia/Shanghai）

## 固定来源

| 资产 | 固定标识 | 用途 | 许可状态 |
|---|---|---|---|
| ReactionT5v2 代码 | `https://github.com/sagawatatsuya/ReactionT5v2.git`，commit `76eb08068e10fe255cae5d563a91e1c1e9abac54` | 官方预处理、推理行为与去重实现参考 | MIT，官方仓库 `LICENSE.txt` |
| ReactionT5v2-forward-USPTO_MIT | `sagawa/ReactionT5v2-forward-USPTO_MIT`，revision `260cf4aa7e6f55b1a9db1e328f48701cd3fd7cfc` | 首日唯一主模型 | 模型卡标记 MIT |
| USPTO_MIT.zip | Google Drive file id `1NEb5DSwqFfg4gm0P2lPb19cqPt96-OSC`；SHA-256 `3FB96BEF12B0E87EF0F0EA243723A861B1FDE161F6CD9DC45FE879EB98891C70` | train/validation 作为参考集合，test 作为评价池 | 归档无许可文件；上游 MolecularTransformer 称其为 open-source，底层 USPTO 政府资料原则上为公共领域。仅用于研究，不再分发数据 |
| USPTO fully mapped split | `wengong-jin/nips17-rexgen@fb7dea369b0721b88cd0133a7d66348d244f65d3` 中 `USPTO/data.zip`；SHA-256 `6D94A136E11F76FE464430CB95D1AE6DB37B6CA352161CA4EDDDD9E6FE76A88A` | 仅用于把同一 480K split 的反应中心/模板信息对齐到主 CSV；对齐失败不得强配 | 仓库 MIT；底层 USPTO 资料按 USPTO 条款署名，不再分发派生数据 |

## 纳排与规范化

- 只纳入反应物与产物均可被 RDKit 解析且主产物明确的记录。
- 保存原始文件、原始行号、原始字符串、规范化字符串和 SHA-256。
- 规范化后完全重复只保留确定性首条；test 与 train/validation 的精确重叠必须为 0，否则来源门失败。
- 主分析保留立体化学；去立体化学仅作为预定义敏感性分析。
- 模板提取失败标记为 `template_failed`，不得当作未见模板。
- 预实验模板固定为 reaction-center 环境签名：由 fully mapped 反应和其中心编辑对推导，主定义半径 1、敏感性半径 0；签名必须去除 atom-map 编号并对编辑和端点排序。它不是标准 RDChiral SMARTS，正式研究需另做模板定义敏感性。

## 暴露等级

- `protected_exact`：仅当官方去重实现及本地审计共同支持精确 test 反应被排除。
- `reference_novel`：相对本合同 train/validation 参考集合新颖，但无法排除 ORD 或分子预训练近邻暴露。
- `unknown_exposure`：证据不足；仅描述，不进入严格 OOD 叙述。

## 抽样冻结规则

- test 清洗后按 train-only 参考模板频率与反应物最近邻分层；validation 保留给后续校准。阈值只使用输入和参考集合确定。
- 固定 2,000 条、随机种子 20260809；目标保证最常见与最罕见联合分箱各至少 300 条。
- 若数据本身无法满足分箱计数，不修改阈值救门槛，记录压力门不可评价或失败。
- 文件哈希、实际列名和最终阈值在首次只读数据审计后追加；在此之前不运行模型。

## 已核验的主数据文件

| 文件 | 字节 | SHA-256 |
|---|---:|---|
| `MIT_separated/train.csv` | 48,191,398 | `4C906E960CD96A0F39E1790314DB642742A0FE8152E4481CC7C81AAC83F4ED14` |
| `MIT_separated/val.csv` | 3,526,860 | `1712F93B057F783FB2FBB4F74A89AB9E302E19230D938C4F74C5AD959DF8FC1D` |
| `MIT_separated/test.csv` | 8,612,805 | `B4D870D0550F6CF10A887DCEE2BCFBAA46DB14DBBDEC734FA0F48B7EC661B5CB` |

## 来源门结果

独立复核确认 test 与 train 原始三元组重叠 739 条、与 validation 重叠 50 条，其中 7 条同时重叠；并集 782/40,000。该结果违反零精确重叠合同，因此本合同在抽样和模型推理前终止，不做去重后续跑。
