# Activity-cliff selective-UQ protocol

本分析合同冻结于 2026-08-09。数据集依据规模、有效率、重复情况和 cliff 元数据选取，选择结果记录在 `dataset_selection.json`。

Activity-cliff 标签复现官方 MoleculeACE 默认实现：Morgan 1024-bit Tanimoto、generic-scaffold fingerprint Tanimoto 或规范 SMILES Levenshtein similarity 任一达到 0.90，并且浓度 fold-change 严格大于 10。预测基线使用的 2048-bit ECFP 不改变 cliff 标签。

分析单位为数据集/靶点。三个 seed 用于构造 ensemble 均值和方差。数据、模型、split、标准化器、距离索引、校准残差和评价路径均可由逐分子清单追溯。

决策依据五项预设指标，数据集、cliff 阈值、分组和 test 评价流程保持固定。

## Strict split construction and balance diagnostic

Before model fitting, each molecule was canonicalized and its valid RDKit single-cut MMP fragmentations were enumerated. The canonical isomeric SMILES of the largest heavy-atom fragment was used as each MMP core key. Molecules sharing a core key or an exact Bemis–Murcko scaffold were connected by union-find, and each connected component was allocated as a unit. Components followed the frozen processed-frame row order, were ordered by the SHA-256 digest of `20260809|<component id>`, and were greedily assigned toward train/calibration/test fractions of 0.6/0.2/0.2. Complete membership and assignments are stored in the split manifest. Canonical SMILES, operational MMP-series links and scaffolds do not cross partitions.

Activity-cliff similarity was inherited unchanged from the pinned code: radius-2, 1024-bit Morgan fingerprints; RDKit `MakeScaffoldGeneric` (falling back to `GetScaffoldForMol`) followed by the same Morgan fingerprint for the generic-scaffold similarity; and canonical-SMILES Levenshtein similarity `1 - distance / max(string lengths)`. The three similarities were unioned at 0.90.

## Frozen-frame provenance and realized allocations

`results/split_provenance.csv` records the SHA-256 digest and row count of each processed frame, together with the digest and 16,528-row count of the externally held `data/split_manifest.csv`. The structure-bearing manifest is not distributed in this repository. The provenance table also reports the original and strict train/calibration/test molecule counts and the largest strict connected component for each dataset. The strict allocations were 2,091/610/616 (CHEMBL214_Ki), 1,713/335/301 (CHEMBL235_EC50), and 1,578/510/510 (CHEMBL236_Ki); the largest strict components contained 1,875, 1,086, and 175 molecules, respectively.

For each matched covariate, the absolute standardized mean difference was calculated as `|mean_cliff - mean_noncliff| / sqrt((var_cliff + var_noncliff)/2)` using finite values and NumPy population variances.
