# P42 assembly-triage v1 归档说明

## 发现的缺陷

首次运行的 `5c_triage_assembly_candidates.py` 对 `supplementary` cluster key 的格式误作三段解析。实际格式为 `supplementary|chr11:start-end:strand`，只有两个 `|` 分段中的第二段含坐标和方向。因此 16 个被标为 `assembly_priority` 的 supplementary rows 在输出中 `anchor_chrom` 与 `anchor_1based` 为空。

## 影响边界

该缺陷不改变 5b 的 51,364 个原始 cluster、fragment/read-end 支持计数，亦不改变 5c 的 tier 计数或 selection policy；它只使 v1 supplementary rows 无法被正确合并进后续 local-assembly window。该 v1 输出没有被用于任何 local assembly、DLD1 comparison 或结论。

## 修复与对照

修正后的现行脚本仍为 `scripts/20261023/5c_triage_assembly_candidates.py`，修正输出写入 `analysis/20261023_fullcohort_repeat_resolve/04_junction_evidence/assembly_triage_v2/`。v1 脚本和全部首次输出原样保留在本目录，供逐项比较。
