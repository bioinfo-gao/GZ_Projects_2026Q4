# P51 研究计划与审阅记录 — Lijian Wu 小鼠 bulk RNA-seq（control vs Crtc2 × 4 处理）

- **Plan Date（创建）**：2026-11-07（文件名 MMDD 永不改；修订记在下方更新记录）
- **项目根目录**：`/home/gao/projects_2026Q4/51_Lijian_Wu_RNAseq`（实体 `/mnt/ex_8T_SSD/projects_2026Q4/51_Lijian_Wu_RNAseq`）
- **授权**：用户 2026-11-07 原话「请按照 …Sample_Registration_Form.xlsx 分析项目51 RNAseq，数据还不全，能分析的先开始，自动完成全过程」⇒ 含「自动」「全过程」，全流程（分析 + 交付目录）授权；**不授权对外发送**。
- **样本表 SSOT**：[`scripts/0_sample_table.csv`](../scripts/0_sample_table.csv)（由 `0_build_sample_table.py` 从送样单 + 磁盘实测生成；本文不复制样本表）

## 更新记录
- 2026-11-07 — 创建；上游 partial 阶段（13/24）已在 tmux `p51_rnaseq` 启动，`p51_watch` 守候到齐后自动跑 full。
- 2026-11-08 01:2x — **本地 Dropbox 客户端停滞**：`dropbox status` 报 "Downloading 78 files"，实测网络接收 ≈0.002 MB/s、目录自 11-07 22:51 起无新文件（27/48 落地，云端 48/48；云端全部文件 `modified_time` = 2026-11-08T00:14Z，目录现显示为独立 namespace，推测该时刻被共享/移动后客户端重协调卡住——**推断，未证实**）。守护进程未动（`kill` 被权限拦，按规则不绕过）。改用连接器一次性直链 + curl 把缺失的 21 个文件（65.9 GiB）下载到 `raw_fastq_direct_dropbox_20261108/`（`scripts/0b_direct_download_dropbox.sh`，tmux `p51_dl`，≈11 MB/s）；`scripts/0c_verify_direct_download.py` 按 Dropbox `content_hash`（4 MiB 分块 SHA-256 再 SHA-256）逐个校验，函数先用本地已同步文件 Control_TNFA_2 R1 做阳性对照、哈希一致。脚本 0 改为每个 mate 优先用 Dropbox 同步副本、否则仅用 hash=ok 的直下副本，新增列 `r1_path`/`r2_path`/`fastq_source`；脚本 1 改读绝对路径。已验证现有 13 样本的 samplesheet 与修改前逐字节相同（`-resume` 缓存不受影响）。**两份副本并存**：Dropbox 后续同步完成后同名文件会出现两份（内容相同，hash 可证），SSOT 里 `fastq_source` 记录实际用的是哪份。
- 2026-11-08 03:14 — 21 个直下 FASTQ 全部 `content_hash=ok`（两次被 Dropbox 链接 ~60 min 断连，`curl -C -` 续传）；24/24。full 04:08–07:07 一次成功（447 cached 复用 / 382 新提交，0 ERROR；日志在 `logs/watch_20261107.log`——watcher 内调用时 `$TMUX` 已设，`2_nextflow.sh` 不走 tee 分支，故无单独 `nextflow_full_*.log`）。09:49 本地 Dropbox 已自行恢复，47/48 落地。
- 2026-11-08 09:3x — **低表达过滤口径变更（§2A filter-impact audit）**：house n−2 规则（≥22/24）在 8 组设计中删掉条件特异基因，含 Cxcl2/Cxcl10/Il6（StepA 后 1,280 个"某组 3 重复均 ≥10"的基因被删，其中 377 个只在 TNF 组表达）。用户裁定「样品数量巨大，不能用 n-2 规则」并给出规则：**本对比两组中 ≥3 个样本 ≥10 counts，或全体 ≥30% 样本 ≥10**。v1(n−2)、v2(某组 3 重复全 ≥10) 原样归档于 `superseded_for_comparison/DE_enrichment/`（README 写明）；v3 = 用户规则，每对比单独检验集合 + `independentFiltering=FALSE` + 集合内 BH，全局并集 12,327 基因。Cxcl1/Ccl2 在全部样本 ≤17 counts，属真实低表达、非过滤误删。
- 2026-11-08 10:1x — **GSEA 排序值变更**：house 排序 sign(LFC)×−log10 p 在 4 个 Crtc2 vs control 对比 GO:BP GSEA ≈0（ORA 数百条）；探索性对照：Crtc2-PGE2 vs control-PGE2 用 Wald stat 得 479 条（min padj 8.7e-08）、top 条目与原排序一致。用户裁定改用 DESeq2 未收缩 Wald `stat`；v3 归档 `superseded_for_comparison/DE_enrichment/v3_..._gsea_signlogp/`；v4 DE 与 v3 逐字节相同（仅结果表加 `stat` 列）。
- 2026-11-08 11:xx — 报告/自审/交付完成：`scripts/6_generate_research_report.R`（9 节；Crtc2 只称 "Crtc2 line"、报告实测 mRNA 0.53–0.67×）、`scripts/8_pre_delivery_audit.R`（149 项；含**独立互证**：用交付 count 矩阵重推每对比检验集合、表内重做 BH——首跑抓出报告措辞未写明"取整后 estimated counts"，已改）、`scripts/9_build_delivery.sh`（无 rm；旧日期报告移入 `superseded_for_comparison/report/`）。交付目录 `custom_research_report_20261108/`（532 文件/132 MB），**内部审阅中、未发送**。描述性"应答保留比例"（PGE2 58%/adenosine 32%/TNF 74%）不是交互检验，报告 §8.3/8.4 写明可按需加做 genotype×treatment interaction。

## 1. 客户送样单要点（原文照录，出处 `docs/Sample_Registration_Form.xlsx`）

| 项 | 送样单内容 | 位置 |
| :--- | :---: | :---: |
| 服务类型 | Section 4 每行 `RNA-Seq (mRNA)`；3A 的勾选框**全部未勾** | A100–A123 / A25–A44 |
| 物种 | Section 4 `Species = Mouse`；3B 勾选框**未勾** | D100–D123 |
| 材料 | `☑ Cell line (established)`；`☑ Normal only`；`☑ Snap-frozen`；Sample Type = `Cell`；0.1 million cells | A55/A63/A65, C/F 列 |
| 细胞系名称/细胞类型 | **未提供** ⇒ 报告写 `Cell type: established mouse cell line (name not specified by client)` | — |
| Crtc2 组的性质 | **未说明**是 KO / KD / 过表达 ⇒ 交付物一律只称 "Crtc2 line"，不写 KO | — |
| 特殊设计 | 3D/3E/Sheet 6 全部未勾 ⇒ 标准两组比较 | A69–A73 |
| 阈值 | 3H：padj 0.05、log2FC 1.0、BH 均写 "Use default" | A87–A90 |
| 生物安全 / 患者来源 | 均未勾 | A94–A95 |

**对比（3G，原文，B 为 baseline）**：
1. crtc2-un vs control-un
2. crtc2-tnf vs control-tnf
3. crtc2-pge2 vs control-pge2
4. crtc2-adenosine vs control-adenosine
5. control-pge2 vs control-un
6. control-adenosine vs control-un

严格只做这 6 个（同客户 proj19 的用户裁定：「严格按客户标注，不额外添加」）。**未列出** control-TNF vs control-un、Crtc2 各处理 vs Crtc2-un —— 不做；若客户后续要，按增量补充轮处理。

## 2. 命名、缩写与映射的判定（假设，须可回溯）

| 项 | 读法 | 判定依据 | 已排除的读法 |
| :--- | :---: | :--- | :--- |
| `PEG2`（Section 4 样本名） | = **PGE2**（prostaglandin E2） | 3G 同一客户写 `pge2`；FASTQ 名 `Control_PGE_*`；与 TNF、adenosine 同属细胞外信号刺激 | PEG = polyethylene glycol（无对应语境） |
| 分析 ID | `control-PGE2-1` 等 | 交付物用规范拼写 + 样本表/报告 Sample Information 并列客户原名 `client_sample_name` | 直接沿用 `PEG2`（会把笔误带进全部交付文件名） |
| FASTQ ↔ 样本 | `Control_N`=control-un-N；`Control_Ad_N`=adenosine；`Control_PGE_N`=PGE2；`Control_TNFA_N`=TNF；Crtc2 同理；重复号一一对应 | 8 组 × 3 名称结构一致；Illumina S 号连续成组（S1–3 un、S4–6 Ad、S7–9 PGE、S10–12 TNFA、S13–24 Crtc2 同序） | 无其它可配方式 |
| `un` | untreated | 与三种处理并列 | — |
| `TNF` | TNF（FASTQ 写 TNFA = TNF-α） | — | — |

⚠ 3H 的 log2FC "Use default (1.0)"：送样单模板默认值是 1.0，而本公司 RNA-seq 标准流程（`/rnaseq` skill、本客户 proj19 交付）用 |log2FC| ≥ 0.263（1.2 倍）。**采用我方标准 0.263**（阈值属我方决定），报告 Methods 写明；并在 DEG 汇总里**同时给出 |log2FC| ≥ 1 的计数**，让客户看到模板默认口径下的数字。

## 3. 数据（截至 2026-11-07 22:49 实测）

- 目录：`/Work_bio/dropbox/Dropbox_Data/Quote_09292601/`（Dropbox 渐进同步中）
- **截至 2026-11-07 22:49 本地为 27 个 FASTQ / 69.9 GiB，13/24 样本 R1+R2 成对**；同时刻 Dropbox 云端该文件夹为 **48 个文件**（24 样本 × R1/R2，list_folder `has_more=false`）⇒ 缺口是同步滞后，不是客户漏传，不需要联系客户。
- 单 lane（flowcell `2533TCLT4` lane 3），PE，每样本单对文件；每样本 GiB 见样本表 `R1_GiB/R2_GiB`；`raw_read_pairs` 在上游跑完后由 `4a_fill_read_pairs.py` 从 FastQC 回填。

## 4. 证据模型 decision table（§2A）

| 研究问题 | 观察单位 | 保留/过滤规则（为什么） | 判定指标 | 允许的结论 |
| :--- | :--- | :--- | :--- | :--- |
| 同一处理下 Crtc2 line 与 control 的转录差异（对比 1–4） | 基因（Salmon gene-level estimated counts） | ① 核糖体 + 非编码/预测基因 regex（排除 Rps6k/Rps19bp1/Rps23rg1 look-alike）—— 防 PCA 被核糖体主导、去掉无 Entrez 的 Gm 基因；② ≥10 counts in ≥22/24 样本 —— 去掉检测不到的基因以免稀释多重检验。**删掉会遗漏**：低表达、细胞类型特异的转录因子（记录过滤前后基因数） | DESeq2 Wald + ashr shrink；padj ≤ 0.05 & \|log2FC\| ≥ 0.263 | 该细胞系、该处理下的基因层面差异；**不能**推断 Crtc2 的 KO/OE 性质（客户未提供），只报告 Crtc2 本身的表达测量 |
| control line 对 PGE2 / adenosine 的应答（对比 5–6） | 同上 | 同上 | 同上 | 该细胞系对处理的转录应答 |
| 数据质量 / 样本身份 | 样本 | 不剔除任何样本，除非有错标证据（PCA 最近邻跨组 + 独立证据） | PCA 组内/组间距离比、最近邻同组率、Crtc2 表达是否与基因型一致 | 是否存在错标/离群；若有，报告里给独立判断 |

- **模型**：24 样本同一 run、同一 lane ⇒ 一个 `~Group`（8 水平）模型，6 个对比都从中提取（与 proj19 相同）。不加交互项（3E 未勾 factorial；客户只要两两比较）。
- **Data universe**：全部 24 样本 / 全部 GENCODE M35 基因进入定量；过滤只发生在 DE 前，数目入 `filter_stats.rds` 并写进报告 §5。

## 5. 流程与脚本（`scripts/`）

| 步骤 | 脚本 | 状态 |
| :--- | :--- | :---: |
| 0 样本表 | `0_build_sample_table.py` | 已跑（每 10 min 由 watcher 刷新） |
| 1 samplesheet + gzip -t | `1_produce_nf-core_Samplesheet.py` | 由 2 调用 |
| 2 nf-core/rnaseq 3.15.1 | `2_nextflow.sh partial\|full`（GRCm39/M35；STAR `--twopassMode None --outFilterMultimapNmax 3`；strandedness auto） | partial 运行中 |
| 3 守候/接续 | `3_wait_and_launch.sh`（tmux `p51_watch`） | 运行中 |
| 4a 读数回填 | `4a_fill_read_pairs.py` | 待 full |
| 4 DESeq2/PCA | `4_run_DE_PCA.R`（Step B = v3 用户规则，见更新记录 11-08） | v3 完成 |
| 5 富集 | `5_run_enrichment.R`（GO/KEGG ORA UP/DOWN/ALL + GSEA GO/KEGG/Hallmark；失败记 ERROR 不静默） | 待 full |
| 6 报告 / 8 自审 / 9 交付 | 读结果后编写 | 待 READY_FOR_REPORT |

- **资源**：独占整机（启动时 load 0.3），executor 56 线程 / 108 GB；STAR 18 cpu × 3。
- **打断安全性**：所有输出写 `output_results/`、`work/`、新日期 `Data_Analysis_YYYYMMDD/`；中途打断可 `-resume`；不触碰任何其它项目。

## 6. 待办 / 风险
- partial 阶段与 full 阶段共享 `work/`，full 用 `-resume` 复用 per-sample 任务；需在 full 日志里核对 cached 任务数 > 0（确认复用真的发生）。
- Crtc2 表达方向需在 DE 后核对（proj19 教训：KO 不一定在 mRNA 层面可见）；结果写入 §7。
- README 提到"初始误建版本保留于 superseded_for_comparison/"——**实测本项目内不存在该目录**（2026-11-07 23:5x 检查），未见其它位置；不影响分析，记录备查。
