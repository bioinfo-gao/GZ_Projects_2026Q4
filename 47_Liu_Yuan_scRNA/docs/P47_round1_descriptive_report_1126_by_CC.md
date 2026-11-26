# P47 Round-1 描述性分析报告（内部，Claude Code）

日期：2026-11-26　｜　作者：Claude Code　｜　物种/细胞：Human（*Homo sapiens*）HEK293T 细胞系，10x Flex（Fixed RNA Profiling），8 样本 / 2 个 GEX pool
性质：**内部审阅稿，非客户交付，不对外发送。** 运行口径见 [P47_round1_descriptive_run_decision_table_1125_by_CC.md](P47_round1_descriptive_run_decision_table_1125_by_CC.md)；全部数字来自 `output_results/round1_descriptive_1125/`（下称 `OUT/`）。

## 0. 结论（先读这三条）

1. **edited vs control 的表达差异在本数据里无法与 pool 批次分开。** 跨 pool 的混杂对比 |log2FC|>1 基因数 = **40**，与 pool 内可识别对比的中位数 **40** 相同（比值 1.0）；同时 Harmony 整合后 pool 混合度仅达随机期望的 0.60。任何 "editing 导致 X" 的说法本轮都没有依据。（`OUT/17_pool_effect_vs_within_pool_effect.csv`、`OUT/14_integration_metrics.csv`）
2. **细胞级检验在本设计下会制造上万个假"显著"基因，已实测证明。** 同条件两样本（N1-Control vs N2-Control）细胞级 Wilcoxon：FDR<0.05 有 **10,850** 个基因；同一样本随机对半对照为 **0** 个。所以检验本身没有失控，这些"显著"来自把细胞当作独立重复（伪重复）。这是本轮不报 p 值的直接证据。（`OUT/18_pseudoreplication_demo.csv`）
3. **技术质量足以支撑后续建模，但有三处需要带着走的混杂**：样本间深度差 2.98×（CBE-Control 最浅，median UMI 26k）；doublet 标记率 15–22%（主分析未剔除，另备敏感性分支）；细胞周期 S 期占比在样本间 55–71%，且与 pool 部分同向。自审 14/14 PASS。

**碱基编辑状态本数据不可判定**（FRP 测的是合成探针，不是原始转录本），本轮未做也不应做。

## 1. 数据与 QC

| 样本 | pool | 公司 cell calling | QC 后细胞 | median UMI | median genes | doublet 标记 % | ambient ρ 中位 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| N1-ABE | pool1_edited | 15,245 | 12,636 | 55,189 | 8,988 | 17.5 | 0.312 |
| N2-ABE | pool1_edited | 19,027 | 15,496 | 62,532 | 9,265 | 21.3 | 0.262 |
| N1-CBE | pool1_edited | 20,706 | 17,339 | 58,247 | 9,235 | 21.2 | 0.295 |
| N2-CBE | pool1_edited | 15,726 | 12,799 | 68,574 | 9,471 | 18.0 | 0.290 |
| N1-Control | pool2_control | 20,780 | 17,565 | 70,004 | 9,319 | 22.3 | 0.353 |
| N2-Control | pool2_control | 20,053 | 17,277 | 67,294 | 9,354 | 22.1 | 0.372 |
| ABE-Control | pool2_control | 12,027 | 9,773 | 52,488 | 8,686 | 15.3 | 0.289 |
| CBE-Control | pool2_control | 12,480 | 10,356 | 34,591 | 7,458 | 14.8 | 0.242 |

- QC 后合计 **113,241** 细胞 × **18,054** 基因（面板 18,082，剔除全零 28 个）。median UMI/genes 为 QC 后值，取自 `OUT/16_pseudobulk_sample_coldata.csv`；cell calling 与 doublet 取自 `OUT/13a_doublet_summary.csv`。
- 过滤用**逐样本 MAD**；每样本每条规则的剔除数见 `OUT/13b_filter_impact_audit.csv`，被剔细胞逐个列于 `OUT/13b_excluded_cells.csv.gz`。
- **doublet 只标记、主分析不剔除**：QC 后细胞中被标记比例为 17–26%（`OUT/13b_two_arms_cellcounts.csv`）。不剔除版和剔除版都保留了细胞数；后续若进入建模，两臂都要跑。
- **MT 基因**：Flex 面板含 12 个 MT 基因（缺 MT-ATP8），其中 11 个在数据中有计数；`pct_counts_mt` 均值 1.11%，各样本 0.72–1.54%（`OUT/12_panel_and_depth_facts.json`）。Flex 的 MT% 不能按常规 3′ 的阈值解读。
- **深度不均衡**：公司 `Reads in cells` 最大/最小 = 2.98×，QC 前 median UMI 比 = 2.36×；CBE-Control 最浅（图 F2）。
- **ambient ρ（decontX）中位数为 0.24–0.37**。Flex 同 pool 内多个样本共享背景，因此 ρ 中包含"来自同 pool 其他样本的转录本"，不能当作单样本 soup 来理解（`OUT/15_ambient_rho_summary.csv`，图 F6）。

## 2. 结构：整合、聚类、细胞周期

- 共 2,000 个 HVG（只做标记，未取子集），Harmony 以 pool 作为批次，Leiden 得到 **11 个簇**，seed 20261125。
- **整合指标**（`OUT/14_integration_metrics.csv`）：跨 pool kNN 占比由 0.147 升到 0.301（若完全混合，期望 0.500），即混合度比值由 0.30 升到 0.60；跨 phase kNN 占比由 0.159 降到 0.157，说明细胞周期结构没有被整合抹掉。**由于 pool 与组别完全共线，残余的 pool 分离不能归因于批次，也不能归因于生物学。**
- **簇的 pool 组成**：簇 3、4、10 中 pool1 占 74–86%，簇 9 中 pool2 占 82%，其余簇接近 50/50（`OUT/14_cluster_by_pool_pct.csv`）。这些簇是"pool 富集簇"，**不是**"编辑诱导的细胞状态"，原因同上。
- **细胞周期**：S 期占比为 55–71%，G2M 为 29–44%（`OUT/14_cellcycle_phase_pct.csv`）。**G1 <1% 是打分规则造成的，不是生物学**：scanpy `score_genes_cell_cycle` 只有在 S 和 G2M 两项得分都 <0 时才判为 G1，所以绝对比例不可解读，只能在样本间做相对比较。S 期占比在 pool2 的 N2-Control 最高（70.8%），在 editor-only 两样本最低（≈55%）。
- 293T 是单一细胞系，**不做细胞类型注释**。

## 3. 描述性效应量（只给 log2FC，不给 p/FDR）

共纳入 12,782 个基因（CPM≥1 且在 ≥4 个样本中满足）；pseudobulk 为 sample-level raw counts 求和（`OUT/16_pseudobulk_sample_counts_genes_x_samples.csv.gz`）。

| 对比 | 范围 | n | \|log2FC\|>0.5 | \|log2FC\|>1 | \|log2FC\|>2 | 中位 \|log2FC\| |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| ABE vs CBE | pool1 内，可识别 | 2 vs 2 | 356 | 40 | 2 | 0.110 |
| guide N1 vs N2 | pool1 内，可识别 | 2 vs 2 | 2 | 1 | 1 | 0.051 |
| editor-only vs gRNA-only | pool2 内，可识别 | 2 vs 2 | 4,208 | 956 | 72 | 0.339 |
| ABE-ctrl vs CBE-ctrl | pool2 内，n=1 | 1 vs 1 | 839 | 42 | 0 | 0.179 |
| edited vs control | **跨 pool，混杂** | 4 vs 4 | 1,059 | 40 | 2 | 0.176 |

来源：`OUT/17_contrast_summary.csv`。逐基因表为 `OUT/17_contrast_*.csv.gz`。

**我的判断（推断，未经检验）：**
- **N1 与 N2 两条 guide 在 pool1 内几乎没有表达差异**（1 个基因 |log2FC|>1；单样本对 baseline 也为 1）。这与"guide 本身不改变全局表达"的设想一致，但每组 n=2 且同 pool，只是描述。
- **pool2 内 editor-only vs gRNA-only 的 956 个基因是本表中最大的信号，但它至少与两件事混在一起**：(a) CBE-Control 的深度只有其他样本的约一半；(b) 这两个 arm 的 S 期占比差 8–16 个百分点（55.1/55.4% 对 63.6/70.8%）。在看到深度匹配（下采样）版本之前，不应把它读成"editor 表达本身的效应"。
- pool2 内 N1-Control vs N2-Control 单样本对有 27 个基因 |log2FC|>1，pool1 内的同类对只有 1 个。说明 pool2 内部的样本间变异大于 pool1，这是后续估计 dispersion 时要考虑的一点。

## 4. 方法学自证：为什么不报 p 值

| 对比 | 细胞数 | FDR<0.05 | FDR<0.05 且 \|lfc\|>0.25 |
| :--- | :---: | :---: | :---: |
| N1-Control 随机对半（校准） | 17,565 | 0 | 0 |
| N1-ABE 随机对半（校准） | 12,636 | 0 | 0 |
| N1-Control vs N2-Control（同条件） | 34,842 | 10,850 | 3,050 |
| N1-ABE vs N2-ABE（同条件） | 28,132 | 9,741 | 1,536 |

随机对半为 0，说明检验没有失控；同条件的两个样本就有上万个"显著"基因，说明细胞级 p 值衡量的是样本间差异，而不是处理效应。本设计每条件 n≤2 且主方向与 pool 共线，因此 sample-level 检验也不成立（decision table §1 的三条理由）。

## 5. 图件（7 个 PDF）

目录：`/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125/figures_1125/`

| 文件 | 内容 |
| :--- | :--- |
| F1_qc_per_sample.pdf | 逐样本 UMI / genes / MT%（MAD 过滤后） |
| F2_depth_imbalance.pdf | 公司 reads/cell 与 median UMI |
| F3a_umap_pool_sample.pdf | UMAP：pool + sample |
| F3b_umap_phase_leiden.pdf | UMAP：细胞周期 + Leiden |
| F4_composition.pdf | 每样本的 phase 与 Leiden 组成 |
| F5_pool_vs_within_pool_effect.pdf | 混杂对比 vs pool 内对比 vs 单样本对 baseline |
| F6_ambient_rho.pdf | 逐样本 decontX ρ |

1126 出图时修正了三处（旧版与旧脚本均已归档到 `superseded_for_comparison/`，旧版不能当作当前版本看）：
- pool 与 sample 原本共用蓝/红，现改为 pool 用灰、样本用 8 个有彩色；F3 拆成 a/b 两张图，避免 leiden 与 sample 在同一张图里撞色。
- G2M 原色 #DDF3DE 在白底上看不见（对比度 1.17），现改为浅灰。
- **F5 旧标题写的是 "Pool/batch effect dominates"，与数据（比值 1.0）矛盾**，现改为由数据生成："Confounded cross-pool contrast is 1.0× the within-pool median — group effect not separable from pool"。

## 6. 下一步（需要客户设计信息才能解锁）

与 Codex 版计划 §5 的问题相同，本轮没有新增：8 个样本的 target/control 定义与独立培养关系；两个 GEX pool 是否为独立 GEM/建库批次、是否有桥接样本；是否另有 DNA/amplicon 数据可用于判定编辑状态。解锁后的第一步是用 `OUT/16_*` 的 pseudobulk counts 建模，并同时跑 doublet 两臂与深度匹配敏感性分析。

## 7. 复核方式

- 自审：`OUT/18_self_audit.csv`（14 条，全部 PASS），日志 `logs/18_audit_1125.log`。
- 确定性：同进程两遍 kNN 逐元素一致；跨进程 kNN 有 0.103% 的元素不同，原因是距离并列时的 tie-breaking；leiden 标签与混合度指标两遍完全一致。
- 出图：`logs/19_figures_1125.log`（STEP19_PASS）；运行链路见 `logs/20_chain_round1_1125.log`。
