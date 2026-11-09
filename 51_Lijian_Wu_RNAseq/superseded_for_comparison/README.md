# P51 superseded_for_comparison

内部归档，**不进客户交付**。每个版本 = 当时的 `Data_Analysis_20261108/` 整目录（原始产物，mv 过来，未改动）
+ 当时的 `scripts/rds_cache/` 副本 + 运行日志 + 脚本快照。当前有效版 = 项目根下的 `Data_Analysis_20261108/`（v4）。

## DE_enrichment/ — 低表达过滤规则（Step B）三个版本

| 版本 | 规则 | 保留基因 | Crtc2-TNF vs control-TNF DEG | 缺陷 |
| :--- | :---: | :---: | :---: | :--- |
| v1_20261108_n-minus-2_low_count_filter | ≥10 counts in ≥ n−2 = 22/24 samples（rnaseq skill house rule） | 10,826 | 5,327 | 8 组×3 设计里，只在 1–2 种处理下表达的基因最多 3–6 个样本达标，永远到不了 22 ⇒ 系统性删掉条件特异基因（Cxcl2/Cxcl10/Il6 等）。用户裁定：「样品数量巨大，不能用 n-2 规则，这相当不合理」 |
| v2_20261108_group_all3reps_filter | ≥1 个组的 3 个重复全部 ≥10 | 12,110 | 5,915 | 无实质缺陷；被 v3 取代（用户给定规则） |
| v3_20261108_pair3_or_30pct_filter_gsea_signlogp | 本对比两组 6 样本中 ≥3 个 ≥10，**或**全体 ≥30%（≥8/24）≥10；每个对比单独限定检验集合、`independentFiltering=FALSE`、集合内 BH | 12,327（全局并集；单对比 11,823–12,172） | 5,905 | DE 无缺陷；**GSEA 排序值** sign(LFC)×−log10 p 尾部过重 ⇒ 4 个 Crtc2 vs control 对比 GO:BP GSEA ≈0（Crtc2-PGE2: 0，min padj 0.063；Wald stat 实测 479，top 条目相同）。用户裁定改用 Wald stat |
| **v4（当前，项目根）** | 同 v3；GSEA 排序值 = DESeq2 未收缩 Wald `stat`（结果表新增 `stat` 列） | 12,327 | 5,905（与 v3 逐字节相同） | — |

依据（v3）：DESeq2 vignette pre-filtering = "≥10 counts in ≥ smallest group size"；Michael Love 对 per-contrast 显式过滤的建议（support.bioconductor.org/p/109669）是 ≥10 reads in ≥3 samples + `independentFiltering=FALSE`；GTEx 用 "≥6 reads in ≥20% samples"。

**脚本快照**：v1 = 当时的原文件（`cp`，未改）。**v2 的 `4_run_DE_PCA.R_RECONSTRUCTED` 是重建，不是当时的文件本身** ——
v2 运行中脚本已被改成 v3，故按会话里 v2 编辑文本反向重建；已核对重建版 Step B 打印格式与 v2 运行日志
`downstream_20261108_0936_4.log` 中 "Step B (group-aware): 26584 -> 12110 | old n-2 rule would keep 10826" 一致。

对照：`diff v1_*/scripts_snapshot/4_run_DE_PCA.R ../scripts/4_run_DE_PCA.R`；
DEG 数并排见各版 `Data_Analysis_20261108/DE_PCA_Results/DEG_summary_all_comparisons.csv`；
`v1_*/explore_group_aware_filter_vs_n-2.csv` 是做决定前在 scratchpad 跑的探索性对比（v2 规则）。

v3 的脚本快照取自 git `HEAD`（baa3970，原文件）；`explore_gsea_rank_metric_test.R` 是排序指标的探索性对照测试。

## docs/ — 研究计划与审阅文档

| 版本 | 内容 | 被取代原因 |
| :--- | :---: | :--- |
| v1_20261108_process_log_with_data_issues/P51_research_plan_and_review_1107.md | 原文件 `mv` 过来，未改动。内容为送样单要点、命名判定（PEG2=PGE2 等缩写的读法与依据）、Dropbox 停滞与直下校验、过滤规则与 GSEA 排序的变更记录、流程状态 | 2026-11-08 用户要求：重写成只讲分析本身的清晰文档，数据问题不进正文。**当前有效版 = `docs/P51_research_plan_and_review_1107.md`（同名，内容已重写）** |

对照：`diff superseded_for_comparison/docs/v1_*/P51_research_plan_and_review_1107.md docs/P51_research_plan_and_review_1107.md`。
