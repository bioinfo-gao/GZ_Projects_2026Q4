# P42 junction-cluster 支持单位：归档版本说明

## 版本性质

本目录保存 2026-10-23 修订前的 `5b_cluster_junction_evidence.py`。该脚本**尚未运行，未产生 cluster 结果**；保留它是为了让支持单位的设计变更可直接比较，而不是删除或覆盖历史代码。

## 已发现的问题

该版本将 `(QNAME, read-end)` 命名为 `molecule_id`，并将它直接累加为 `independent_molecule_count`。在 paired-end WGS 中，同一 QNAME 的 R1 和 R2 来自同一个 library template/fragment；若两端均落入同一 cluster，该口径会把一个 fragment 放大为两个独立支持单位。

## 替代版本与用户确认

用户于 2026-10-23 确认：以同一样本内唯一 `(QNAME/template)` 作为独立 fragment 的主阈值，`QNAME+read-end` 只作为并列 read-level 支持数。现行脚本位于 `scripts/20261023/5b_cluster_junction_evidence.py`，输出 `independent_fragment_count` 与 `independent_read_end_count` 两列，并保留所有 raw record context。

## 对照范围

本归档仅涉及计数单位定义；不改变 discovery window、输入 evidence records、cluster binning 规则或候选并非 SV call 的解释边界。
