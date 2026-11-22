# P47 开跑前缺陷清单

Plan Date（创建）: 2026-11-22
状态：内部文档。仅盘点缺陷，未启动正式分析，未发送任何邮件。

## 0. 结论

**环境基本就绪，但还不能按最佳实践完整分析 P47。** 阻塞分四类：客户输入未到（§1）、最佳实践在 P47 上有做不满的地方（§2）、还缺两个工具（§3）、开跑前有几项口径要用户确认（§4）。其中最关键的是 §1：靶点在 mtDNA 还是核基因、8 套样本是否独立、交付的是 FASTQ/BAM 还是只有矩阵，这三点决定后面每一步能做到哪一层。

建议顺序（§5）：工具已补齐（1122 装好 ClusterDE、countsplit；xgboost 实测定为 3.3.1.1）→ 数据到后先做靶位点覆盖审计 → 按覆盖结果定稿 decision table → 用户确认后再启动。

## 1. 客户输入未到

来源：[计划](P47_research_plan_and_review_1029.md) §6；[调研](P47_scRNAseq_best_practices_github_survey_1104.md) §8 末段。样本表 [../scripts/0_sample_table.csv](../scripts/0_sample_table.csv) 截至 2026-11-22 只有表头。数据预计 2026-11-24（周二）到达，尚未定位/清点。

| 缺的输入 | 卡住什么 |
| :--- | :---: |
| 靶点是 mtDNA 还是核编码基因；坐标、参考碱基、预期转换 | 等位基因分支的工具选择与可达深度 |
| 8 套样本是独立培养、克隆还是技术重复；batch | n 的计算、能否做跨条件推断 |
| base editor control 与 gRNA control 各含什么构建体、应怎么比较 | 比较关系与设计矩阵 |
| 交付形式：FASTQ / BAM（含 CB/UB 标签）/ 仅表达矩阵 | 只有矩阵 ⇒ 碱基层面证据无法做，只能报告"不可评估" |

## 2. 最佳实践在 P47 上做不满的地方（设计/数据固有，不是环境问题）

来源：调研 §5.3、§6。

| 环节 | 限制 | 处理 |
| :--- | :---: | :---: |
| 跨条件统计 | 两类对照各 2 套，低于每组 n≥3 门槛 | 只报描述性效应量与一致性，不出正式显著性结论 |
| doublet 检测 | 293T 为同质细胞系；本机阳性对照中 scDblFinder 同质群 recall ≈0.41–0.46（异质群 0.72–0.93） | 仍运行，结合文库大小上尾报告，Limitations 写明漏检 |
| mtDNA 靶点覆盖（若靶点在 mtDNA） | **1122 在公开 10x 3′ v3 pbmc_1k_v3（~54k reads/细胞，血细胞，非 293T）实测**（chrM 每 100 bp 166 个位点，分子层面）：每位点有 ≥1 分子的细胞比例中位 70%（IQR 41–94%）；有覆盖的细胞在该位点中位只有 2 个分子（IQR 1–4.75）；≥90% 细胞有覆盖的位点仅 50/166。⇒ 覆盖高度依赖位点，多数细胞每位点 1–2 分子。[../pipeline_validation_1122/mt_feasibility/mt_coverage_summary_1122.log](../pipeline_validation_1122/mt_feasibility/mt_coverage_summary_1122.log)。293T 的 mtDNA 转录量与测序深度不同，只能作量级参考 | 先做样本级 allele 比例，覆盖审计后再判断能否下沉到单细胞；无覆盖记 NA，不记"未编辑" |
| ambient RNA 校正 | 单一培养物 ambient 谱与细胞谱接近，可辨识性低 | 默认不校正进 DE，只估计污染水平并做敏感性分析 |
| 亚群声称 | 同质细胞系易被切出假亚群 | 声称亚群前须做校准（ClusterDE 或 count splitting，见 §3） |

## 3. 环境与工具缺口

**已就绪（2026-11-21 实测）**：

| 环境 | 内容 | 验证程度 |
| :--- | :---: | :---: |
| `/mnt/ex_8T_SSD/conda_envs/sc_tut_main` | scanpy 1.11.5（内置 `sc.pp.scrublet`）、scvi 1.3.3、harmonypy、pertpy、decoupler | 教程已验证 |
| `/mnt/ex_8T_SSD/conda_envs/p47_r` | SoupX 1.6.2、celda 1.26.0、scDblFinder 1.24.10 + xgboost 3.3.1.1、DESeq2、edgeR、zellkonverter、DropletUtils；**1122 增补** countsplit 4.0.0、scDesign3 1.7.6、ClusterDE 0.99.4（Seurat 5.5.1 已作为依赖在内） | 加载 + 玩具数据冒烟通过：decontX 盲估污染 0.093（真值 0.10），[check6.log](P47_scRNAseq_best_practices_evidence_1104/check6.log)；ClusterDE/countsplit 阴性对照（同质群硬切 2 簇）naive 64 个假阳性 → 两者均 0，阳性对照真 DE 命中 countsplit 20/20、ClusterDE 19/20，[check5.log](P47_scRNAseq_best_practices_evidence_1104/check5.log) |
| `/mnt/ex_8T_SSD/conda_envs/p47_allele` | cellsnp-lite 1.2.3、maegatk 0.2.0、samtools 1.24、pysam 0.24.1 | cellsnp-lite 合成 BAM 阳性对照通过；**maegatk 只验证到加载** |

建法与校验：[../scripts/1_setup_envs.sh](../scripts/1_setup_envs.sh)（`--verify-only` 可复核）；日志在 `../logs/setup_envs_*.log`。

**仍缺**：

| 缺口 | 实测 | 影响 | 建议 |
| :--- | :---: | :---: | :---: |
| ~~ClusterDE、countsplit~~ | **1122 已装进 `p47_r` 并过阴/阳性对照**（见上表） | — | 已解决；用法坑见下 |
| ~~Cell Ranger~~ | **1122 已用 STARsolo 替代并验证**：STAR 2.7.11b（`/mnt/ex_8T_SSD/conda_envs/p47_align`）+ [../scripts/4_starsolo.sh](../scripts/4_starsolo.sh)；在 10x 公开 pbmc_1k_v3 上与 10x 官方 Cell Ranger 3.0 结果逐细胞对比：细胞重叠 99.8%、每细胞 UMI Spearman 0.9998、每细胞表达谱 Pearson 中位 0.971（[../pipeline_validation_1122/compare_starsolo_vs_cellranger_1122.log](../pipeline_validation_1122/compare_starsolo_vs_cellranger_1122.log)） | — | 已解决；不需要装 Cell Ranger |
| ~~10x 参考包 / 白名单~~ | **1122 已建**：`/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/starsolo_10x2024Afilter_star2711b/`（GENCODE v45 + 按 10x 2024-A 规则过滤，39,116 基因）；白名单 `/Work_bio/references/10x_whitelists_ngs_tools_1.8.6/`（v3 在 pbmc_1k_v3 上精确匹配 95.0%） | — | 已解决；化学由 [../scripts/3_detect_chemistry.py](../scripts/3_detect_chemistry.py) 从数据判定（客户邮件未提化学） |
| maegatk 功能适配 | 0.2.0 实测：`--min-reads` 默认 **1**（调研说的"<3 丢弃"与本版本不符）、`--NHmax` 2；**不加 `-so` 就失败但退出码 0**（dash 解析 `&>` 的 bug，详见 `~/Utilities/software_pitfalls/maegatk.md`）；pbmc chrM 全量功能测试见 [../pipeline_validation_1122/mt_feasibility/](../pipeline_validation_1122/mt_feasibility/) | 用错调用方式会静默得到空结果；**depth 不是分子数** | 永远 `-so` + 先 `mamba activate`。**1122 全量跑通**（1,226/1,227 细胞），与审计脚本位点间 Spearman 0.994、位点×细胞 0.956，但 maegatk depth 是 (CB,UB) 分子数的中位 1.53 倍：它按 (CB,UB,read 起点) 计数（位点 12850：27,508 vs maegatk 27,450；(CB,UB) 仅 11,933）。⇒ **主计数用 [../scripts/5_target_site_audit.py](../scripts/5_target_site_audit.py) 的 (CB,UB) 分子；maegatk 只作算法对照，不用其 depth 当 n**。[../pipeline_validation_1122/mt_feasibility/compare_maegatk_vs_audit_1122.log](../pipeline_validation_1122/mt_feasibility/compare_maegatk_vs_audit_1122.log) |
| scib-metrics | `sc_tut_main` 无 | 只影响整合效果量化；P47 主线不整合 | 暂不装 |

**ClusterDE / countsplit 用法坑（1122 冒烟时踩到，均不报错或只给警告）**：

| 坑 | 后果 | 正确做法 |
| :--- | :---: | :---: |
| `constructNull()` 只接受 **Seurat 对象**（源码直接调 `Seurat::GetAssayData`），Seurat 却只在 Suggests | 传 SCE 报 `UseMethod GetAssayData` | 先 `CreateSeuratObject()` |
| 0.99.4 的 `callDE()` 返回 tibble（`gene` / `cs` / `record`），没有 `$DEgenes` | 取 `res$DEgenes` 只给一条 warning，**静默得到 0 个 DE** | `res$gene[res$record == 1]`；断言 `nrow(res) == 基因数` |
| `countsplit()` 要求 **cell × gene** | 不给 `overdisps` 时 Poisson 拆分逐元素进行，方向反了也不报错 | 先转置；给 `overdisps` 时它会检查列数 |
| 不给 `overdisps` ⇒ 默认 Poisson 拆分 | 对负二项数据 train/test 不独立，DE 虚高（本机阳性对照 44 vs 给定后 33 个调用） | 先估每基因 overdispersion 再传入 |

## 4. 开跑前需用户确认的口径

| 事项 | 选项 | 我的倾向与依据 |
| :--- | :---: | :---: |
| ~~scDblFinder 的 xgboost 版本~~ **已定：3.3.1.1（1122 实测裁定）** | — | **更正：此前"倾向 3.4.2.1"是错的。** 直接调用 `scDblFinder:::.xgbtrain`：3.4.2.1 报 `'early_stopping_rounds' requires passing 'eval_set'`，被 `.scDblscore` 的 `tryCatch(..., error = function(e) d$score)` **静默吞掉、退回训练前的初始分数**；3.3.1.1 训练正常。⇒ 1104 那次 3.4.2.1 的较高 recall 不是 xgboost 分数，而是回退分数，两版 recall 对比不能用来选版本。证据：[check7_scdblfinder_xgbtrain_direct.R](P47_scRNAseq_best_practices_evidence_1104/check7_scdblfinder_xgbtrain_direct.R)、`check7_p47_r.log` / `check7_sc_tut_r.log`。正式运行前每次先跑 check7 断言 TRAIN_OK。原对比数据保留：3.4.2.1：[check4.log](P47_scRNAseq_best_practices_evidence_1104/check4.log)（1104，`sc_tut_r`）；3.3.1.1：[check4_p47r_xgboost331_1121.log](P47_scRNAseq_best_practices_evidence_1104/check4_p47r_xgboost331_1121.log)（1121，`p47_r`，同一脚本） |
| 上游定量 | ~~STARsolo / Cell Ranger~~ **STARsolo（1122 验证后定）** | 客户交 FASTQ ⇒ 本机 STARsolo；客户交 Cell Ranger 输出且含 BAM（CB/UB）⇒ 可直接用其 BAM，但建议仍用 STARsolo 重跑一份做 mapping 可重复性对照（只证算法层面） |
| **新增：NUMT 区的 MAPQ 口径**（§2A 要求开跑前确认） | 仅 MAPQ≥20 primary / 另报 NH=2 的 read | **实测（pbmc chrM 166 个网格位点）**：部分区域 45–90% 的 record 是 secondary（chrM 5950、4550–4650、7050–7450 等），同处 primary 中又有一大半因 MAPQ<20 被滤（STAR 给 NH=2 的 read MAPQ 3）。⇒ 若靶点落在这类区域，"MAPQ≥20"会丢掉一半以上证据。**建议**：主表用 MAPQ≥20 primary；同时报告 NH=2 且另一处落在核 NUMT 的分子作为单独一层，不并入主表。客户给出靶点坐标后先查该位点 secondary 比例再定 |
| 计划 §3 decision table 定稿 | — | 须在覆盖审计（§5 第 2 步）后定，按规则 §2A 正式运行前确认 |
| 调研 §8 未执行的三项 | 更新 `/scrna` skill 4 处方法学描述；修订计划 §2 maegatk 行；scanpy 是否升级 | skill 与计划修订建议做；scanpy 固定 1.11.5 写进 Methods |

## 5. 建议顺序

1. **现在**：~~ClusterDE、countsplit 装进 `p47_r`~~（1122 已完成）；~~把 xgboost 版本定下来~~（1122 实测定为 3.3.1.1，见 §4）。
2. **数据到后**：按"截至 YYYY-MM-DD HH:MM 目录内 N 个文件"清点；先做靶位点覆盖与 barcode/UMI 标签审计（计划 §4 第 2 步），**不先跑完整表达流程**。
3. **按覆盖结果**定稿 decision table，用户确认后，以"全流程"启动正式分析。
