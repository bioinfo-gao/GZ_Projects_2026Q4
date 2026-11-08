# P51 superseded_for_comparison

内部归档，**不进客户交付**。每个版本 = 当时的 `Data_Analysis_20261108/` 整目录（原始产物，mv 过来，未改动）
+ 当时的 `scripts/rds_cache/` 副本 + 运行日志 + 脚本快照。当前有效版 = 项目根下的 `Data_Analysis_20261108/`（v3）。

## DE_enrichment/ — 低表达过滤规则（Step B）三个版本

| 版本 | 规则 | 保留基因 | Crtc2-TNF vs control-TNF DEG | 缺陷 |
| :--- | :---: | :---: | :---: | :--- |
| v1_20261108_n-minus-2_low_count_filter | ≥10 counts in ≥ n−2 = 22/24 samples（rnaseq skill house rule） | 10,826 | 5,327 | 8 组×3 设计里，只在 1–2 种处理下表达的基因最多 3–6 个样本达标，永远到不了 22 ⇒ 系统性删掉条件特异基因（Cxcl2/Cxcl10/Il6 等）。用户裁定：「样品数量巨大，不能用 n-2 规则，这相当不合理」 |
| v2_20261108_group_all3reps_filter | ≥1 个组的 3 个重复全部 ≥10 | 12,110 | 5,915 | 无实质缺陷；被 v3 取代（用户给定规则） |
| **v3（当前，项目根）** | 本对比两组 6 样本中 ≥3 个 ≥10，**或**全体 ≥30%（≥8/24）≥10；每个对比单独限定检验集合、`independentFiltering=FALSE`、集合内 BH | 12,327（全局并集；单对比 11,823–12,172） | 5,905 | — |

依据（v3）：DESeq2 vignette pre-filtering = "≥10 counts in ≥ smallest group size"；Michael Love 对 per-contrast 显式过滤的建议（support.bioconductor.org/p/109669）是 ≥10 reads in ≥3 samples + `independentFiltering=FALSE`；GTEx 用 "≥6 reads in ≥20% samples"。

**脚本快照**：v1 = 当时的原文件（`cp`，未改）。**v2 的 `4_run_DE_PCA.R_RECONSTRUCTED` 是重建，不是当时的文件本身** ——
v2 运行中脚本已被改成 v3，故按会话里 v2 编辑文本反向重建；已核对重建版 Step B 打印格式与 v2 运行日志
`downstream_20261108_0936_4.log` 中 "Step B (group-aware): 26584 -> 12110 | old n-2 rule would keep 10826" 一致。

对照：`diff v1_*/scripts_snapshot/4_run_DE_PCA.R ../scripts/4_run_DE_PCA.R`；
DEG 数并排见各版 `Data_Analysis_20261108/DE_PCA_Results/DEG_summary_all_comparisons.csv`；
`v1_*/explore_group_aware_filter_vs_n-2.csv` 是做决定前在 scratchpad 跑的探索性对比（v2 规则）。
