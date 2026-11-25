# P47 — Liu Yuan：293T 单细胞编辑效应分析（研究计划与审阅，Claude Code 版）

Plan Date（创建）: 2026-11-25　｜　作者：Claude Code　｜　性质：内部计划文档，非客户交付物

本文是对同目录 [P47_research_plan_and_review_1125.md](P47_research_plan_and_review_1125.md)（Codex 撰写）的复核与增补版本。**两份并存，Codex 版未被取代、其主体判断经复核成立**；本文记录 4 项新实测证据、2 项实质修正和 1 项安全提示。差异清单见 §7。

## 0. 结论（三条，均为实测）

1. **数据完整，文件数正确。** 本地 380/380 文件、FASTQ md5 8/8 OK。"8 个样却只有 4 对 R1/R2" 不是缺文件：这是 10x **Fixed RNA Profiling（Flex）** 的多重化设计，8 样 = **2 个 GEX library × 每 library 内 4 个 probe barcode**，样本标签藏在 **R2 读段内部**而不在文件名或 index read 里。已从原始 reads 直接验证 8 个 probe barcode 全部在位、~96% reads 可归属。完整格式实测见 [data_format/P47_10x_Flex_data_format_1125.md](../data_format/P47_10x_Flex_data_format_1125.md)。

2. **碱基编辑判定在本批数据上不可行，且阻断不止一重。** 实测 R2 前 50 nt 有 **88.5% 与合成探针序列逐字精确匹配**（正向，反向互补 0 命中）—— 测序读到的是探针 DNA，不是转录本碱基。**新增的更严限制：** 探针面板对 chrM 只覆盖 **600 bp / 3.62%**，**mt-tRNA 与 mt-rRNA 探针数为 0**，且 MT 基因**每基因仅 1 条探针**（无基因内对照）。六个经典致病位点（m.3243/3460/8344/8993/11778/14484）**无一落在任何探针足迹内**。⇒ 即便是"探针信号下降"这条间接路线也基本不可用。**索取 BAM 不能解除此限制**（BAM 里比对的也是探针）。

3. **表达分析不是全盘不可做，但可做的对比比原先设想的窄。** edited 与 control **完全落在不同 pool**，该方向的组间差异与 pool/批次效应完全共线、不可识别 —— 此点 Codex 版判断正确。**但 pool 内部的对比不受此混杂影响**（如 pool 1 内 ABE vs CBE），这类对比是本轮唯一可识别的差异分析入口。另需处理两个新发现的混杂：样本间 read 分配差 **2.98 倍**，且复杂度最低的样本（CBE-Control，26,250 UMI/cell）正是对照组成员。

**当前状态：正式下游分析未启动，这不是阴性结果。** 阻塞项是样本设计定义（§5），不是技术能力。不对外发送任何材料。

## 1. 客户目标、数据状态与来源

客户来信（[client_materials/1_initial_introduction.txt](../client_materials/1_initial_introduction.txt)，Liu Yuan，2026-10-28）原文要点：8 套 293T scRNA-seq，含 **4 个 mitochondrial gene-edited + 2 个 base editor control + 2 个 guide RNA control**；主要目标为 *"determine whether the gene-editing effects can be detected in the single-cell RNA-seq data"*。

"mitochondrial gene-edited" 尚不能确定指 mtDNA 编辑还是核编码线粒体相关基因。**这个区分现在变得关键**：若为 mtDNA 编辑，§2.2 的面板覆盖率问题直接决定靶基因是否被定量；若为核编码基因，覆盖通常不成问题（核基因多为 3 探针/基因）。

### 1.1 下载与完整性（复核通过）

数据位于用户指定目录：

`/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/genewiz_30-1371703644_1124/30-1371703644/`

| 检查项 | 远端清单 | 本地实测 | 判定 |
| :--- | :---: | :---: | :---: |
| 全交付目录文件数 | 380 | 380 | 一致 |
| `00_fastq/` FASTQ 数 | 8 | 8 | 一致 |
| FASTQ md5 | — | 8 OK / 0 FAILED | 通过 |
| `00_fastq/` 体积 | 882.3 GiB | 882.345 GiB | 一致 |
| 全目录体积 | 996.7 GiB | ~997 G | 一致 |

- 远端清单：`/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/remote_listing_genewiz_20261124_0017.txt`（650 行），其 `00_fastq/` 段落只列出这 8 个 FASTQ + 8 个 `.md5`。
- 看门狗日志 `/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/logs/download_watch.log` 末行：`[2026-11-25 01:47:18] DONE 下载结束。本地文件 380（远端清单 380），md5 OK 8/8，FAILED 0。自动重启 0 次。`
- 记录口径：**截至 2026-11-25 01:47:18，目录内为 380 个文件**（Dropbox 式渐进同步不适用于此处，本批为 sFTP 一次性拉取且远端清单已比对一致）。
- 小更正：Codex 版写 `00_fastq/` 有 16 个文件 —— 数量正确，但 8 个 `.md5` 实际位于 `00_fastq/md5/` 子目录，`00_fastq/` 根下只有 8 个 FASTQ。这是用户看目录时觉得"奇怪"的次要成因之一。

### 1.2 公司分析流程与 8 样本拆分（已从原始 reads 独立验证）

两份 `config.csv` 一致记录：Cell Ranger **8.0.0**、参考 `refdata-gex-GRCh38-2020-A`、探针面板 `Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv`，按 **Probe Barcode** 拆样。本地镜像**不含任何 `.bam`/`.cram`**（尽管 `config.csv` 写有 `create-bam,true`）。

**8 样本↔2 文件的对应关系（probe-BC 序列由本地 reads 实测反查，非公司声明）：**

| 样本 ID | 暂定组别 | GEX pool | Probe BC ID | probe-BC 序列 | 该 pool 前 200k reads 占比 | GEX Cells |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| N1-ABE | edited | pool 1 | A-A10 | GAAGGTCCCA | 19.0% | 15,245 |
| N2-ABE | edited | pool 1 | A-B10 | TGCAGTAACT | 25.7% | 19,027 |
| N1-CBE | edited | pool 1 | A-C10 | AGTATGGCTG | 27.6% | 20,706 |
| N2-CBE | edited | pool 1 | A-D10 | ACGACCCATC | 23.8% | 15,726 |
| N1-Control | control | pool 2 | A-E10 | TGATCATGCA | 34.8% | 20,780 |
| N2-Control | control | pool 2 | A-F10 | GAGGGCAAAT | 34.6% | 20,053 |
| CBE-Control | control | pool 2 | A-G10 | ACTACTCGAG | 11.8% | 12,480 |
| ABE-Control | control | pool 2 | A-H10 | AGCTGGCTAC | 15.0% | 12,027 |

`cells_per_tag.json` 列全 384 个 probe barcode 位，每 library 内**仅 4 个非空**，细胞数与 QC 表逐一吻合。合计 136,044 cells。复现：`conda run -n regular_bioinfo python3 -I /home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/data_format/evidence_scripts/2_probe_barcode_demux_check.py`

样本级 QC 指标与 FASTQ/pool 级指标见 [data_format 文档 §5](../data_format/P47_10x_Flex_data_format_1125.md)，不在此重复（SSOT）。

## 2. 可行性与证据边界

### 2.1 当前 FRP 数据不能判定编辑碱基（Codex 版结论成立，此处换成本地实测证据）

Codex 版引 10x 官方文档作为依据，结论正确。本文改以**本地实测**作为依据，证据强度更高：

取 pool 1 的 R2 前 20,000 条 read，将 bases 1–50 与面板 54,580 条探针序列做精确比对：

| 比对方式 | 命中 | 占比 |
| :--- | :---: | :---: |
| 正向精确匹配探针序列 | 17,692 | **88.5%** |
| 反向互补匹配 | 0 | 0% |
| 未精确匹配（多为含测序错误的同批探针） | 2,308 | 11.5% |

面板内 54,580 条探针**长度均为 50 nt**；R2 实测 90 nt = 50 nt 探针 + 常量 capture/linker + 10 nt probe barcode。

**⇒ 测序信息来源是合成 DNA 探针。原转录本上的碱基变异不进入数据。** 本批数据在原理上不能产出编辑位点 ref/alt、编辑率/heteroplasmy，或 base-level 的"已编辑/未编辑"细胞分型。**此结论与测序深度无关，也不能通过索取 BAM 解除。**

### 2.2 【新增】探针面板对 chrM 的覆盖是第二重、且更硬的阻断

这是 Codex 版未涉及、但实质改变可行性判断的一层。实测面板 `probe_set.csv`（54,580 探针 / 19,074 基因）：

| 项目 | 实测 |
| :--- | :---: |
| 面板内 `MT-` 基因数 | **12**（全部蛋白编码） |
| mt-tRNA（`MT-T*`）探针 | **0** |
| mt-rRNA（`MT-RNR1/2`）探针 | **0** |
| 每个 MT 基因探针数 | **1**（全部 `unspliced`） |
| `included=FALSE` 的 MT 探针 | MT-CO1 |
| chrM 被探针覆盖总量 | **600 bp = 16,569 bp 的 3.62%**（仅 included：550 bp = 3.32%） |

12 条探针全部唯一定位到 chrM，足迹与连接点坐标见机器可读 SSOT [data_format/chrM_probe_footprints.csv](../data_format/chrM_probe_footprints.csv)（12 行）。

**六个经典致病 mtDNA 位点逐一检查，无一落在任何探针足迹内：**

| 位点 | 常见关联 | 是否在探针足迹内 |
| :--- | :---: | :---: |
| m.3243A>G | MT-TL1（mt-tRNA），MELAS | **否** |
| m.3460G>A | MT-ND1，LHON | **否** |
| m.8344A>G | MT-TK（mt-tRNA），MERRF | **否** |
| m.8993T>G | MT-ATP6，Leigh | **否** |
| m.11778G>A | MT-ND4，LHON | **否** |
| m.14484T>C | MT-ND6，LHON | **否** |

**三重阻断汇总：**

| 阻断 | 依据 | 后果 |
| :--- | :---: | :---: |
| ① 落点概率低 | 覆盖 3.62% | 编辑位点先验上 ~96% 概率不在任何探针下 |
| ② 经典靶点零覆盖 | 上表 6/6 未覆盖 | 若靶点为 mt-tRNA，该基因**根本不被定量** |
| ③ 无基因内对照 | MT 基因 1 探针/基因 | "探针信号降低"与"该基因表达真降低"完全共线，不可分离 |

第 ③ 条是方法学上最要命的：核基因多有 3 条探针，可用"同基因其他探针作内参"的 probe-ratio 把杂交失败从表达变化中分离；**MT 基因没有这个内参**。因此**即使编辑位点恰落在足迹内，单探针下降也不能定性为编辑证据**，只能作为 hypothesis-generating 的探索性观察。

**⇒ 操作结论：拿到客户靶点坐标后，第一个动作是查它是否落在那 600 bp 内。** 一条命令即可判定（见 [data_format 文档 §3.3](../data_format/P47_10x_Flex_data_format_1125.md)）。若不在，这条间接路线连探索都无从开始；若在，仍只能作探索性呈现。

### 2.3 表达组间比较：edited vs control 方向确实不可识别，但 pool 内对比可用

**实测事实：** pool 1 = 4 个 edited（A-A10…A-D10），pool 2 = 4 个 control（A-E10…A-H10）。两 pool 为不同 physical GEX library、不同 FASTQ、不同测序 run 统计。**无任何样本跨 pool 出现，故无桥接样本。**

**推断（我的判断，需用户/客户确认）：** 两 pool 对应不同 GEM 与建库批次。若成立，则：

| 对比方向 | 与 pool 的关系 | 可识别性 |
| :--- | :---: | :---: |
| edited vs control（4 vs 4） | **与 pool 完全共线** | **不可识别**，统计校正无法从完全共线设计中分离处理效应 |
| pool 1 内：ABE vs CBE（2 vs 2） | pool 内，无 pool 混杂 | **可识别** |
| pool 1 内：N1 vs N2（2 vs 2） | pool 内 | **可识别**（若 N1/N2 确为不同靶点/gRNA） |
| pool 2 内：各 control 两两 | pool 内 | **可识别** |
| 跨 pool 的任意单样本配对 | 与 pool 共线 | **不可识别** |

**这是对 Codex 版的实质修正方向：** Codex 版把表达差异分析整体判为"当前不能归因"，就 edited-vs-control 主方向而言正确，但**结论过宽** —— pool 内部的对比不受该混杂影响，是本轮唯一可识别的差异分析入口，应保留而非一并写掉。

### 2.4 【新增】两个必须在分析期处理的混杂

**(a) 样本间 read 分配差 2.98 倍。** `Reads in cells` 从 7.47e8（CBE-Control）到 2.22e9（N1-Control）。这是 probe barcode 间上样量/捕获效率差异，不是生物学差异。⇒ 样本间比较 UMI/genes per cell 之前必须深度校正或下采样。

**(b) 复杂度最低的样本是对照组成员。** CBE-Control 实测 26,250 UMI/cell、6,840 genes/cell，均为 8 样本最低（其余样本 42,071–61,946 UMI）。⇒ **若后续出现"CBE 相关对比差异最大"的结果，必须先排除这是深度/质量 artifact 而非生物学效应。** 这是本项目最可能产生假阳性结论的单点。

**(c) 同 pool 内样本间存在结构性交叉污染。** Flex 为固定后混样再上机，同一 GEM 液滴可含多个样本的细胞。实测 `frp_gem_barcode_overlap.csv`：同 pool 内任意两 probe barcode 的 GEM 重叠率 **11.4%–17.4%**。⇒ 同 pool 样本并非完全独立观测；ambient 背景校正须按"pool 内多样本"情形处理，不能当单样本 soup 处理。

### 2.5 不把常规 STARsolo 流程用于 Flex（Codex 版结论成立，补充结构性理由）

项目 `scripts/` 下 `2_build_starsolo_ref.sh`、`3_detect_chemistry.py`、`4_starsolo.sh` 针对常规 10x v2/v3/v4 UMI 文库编写，**对 Flex 不适用，且不是调参能解决的**：

| 常规流程假设 | Flex 实测 |
| :--- | :---: |
| cell barcode 16 nt、全在 R1 | **26 nt**，其中 10 nt 在 **R2 内部** |
| R2 是 cDNA，比对到基因组/转录组 | R2 是 50 nt 合成探针，需比对到 **probe set** |
| 用 10x whitelist 校正 barcode | 需同时校正 GEM whitelist **和** 384 个 probe barcode |
| 计数单位 = 基因 | 计数单位 = **探针**，再汇总到基因 |

建议把这三个脚本**标注为"仅适用于常规 3′/5′ 文库，本项目不使用"而非删除**（它们在 1122 的管线验证中已跑通，是可复用资产）。

## 3. Decision table：本轮数据可支持的分析

| 研究问题 / data universe | 观察或证据单位 | 保留、去重与过滤规则 | 判定指标与 filter-impact audit | 当前允许的结论 |
| :--- | :---: | :---: | :---: | :---: |
| 编辑碱基是否出现；全部 FRP FASTQ 与公司产物 | 需原始转录本碱基或 DNA molecule；实际读到的是 synthetic probe | 不把 probe sequence / molecule-info 解释为 transcript sequence；不按表达阈值定义"已编辑细胞" | 标为 assay-level **not assessable**；不生成 ref/alt、编辑率或"未编辑"计数 | 本批 FRP 不可用于碱基级编辑判定；**结论不是"没有检测到编辑"** |
| 靶点是否落在探针足迹内（前置可行性闸门） | chrM 1-based 坐标 vs 12 条探针足迹 | 以 `chrM_probe_footprints.csv` 为唯一判据；`included=FALSE` 的 MT-CO1 单列 | 输出命中/未命中及所属基因、连接点距离 | 可给出二分判定；未命中 ⇒ 间接路线不成立，分析不启动 |
| 每样本表达与文库 QC；8 份公司 raw/filtered 矩阵 | sample × cell；跨条件单位为独立培养样本 | filtered = 公司 cell-called 集合；raw 用于 cell-calling/ambient 核查；逐样本记录 barcode 保留/排除原因 | 每样本 cells、UMI/基因分布、mapping/probe metrics；输出过滤前后 barcode 数及 excluded table | 可报告供应商 QC 与本项目描述性表达 QC；**必须同时报告 read 深度差 2.98×** |
| pool 内可识别对比（ABE vs CBE 等） | 独立生物学样本 × cell state × gene 的 raw-count pseudobulk；细胞数不充当重复数 | 仅在同 pool 内构造 contrast；不合并不同 control 类型；不用 cell-level p 值伪增样本量 | 报告 n、设计矩阵秩/共线性、过滤影响；对 CBE-Control 低复杂度做敏感性分析 | **可做探索性 pseudobulk 对比**，n 极小须显式声明 |
| edited vs control 的表达响应 | 同上 | 同上 | 先审 biological replicate 与 pool 交叉；**实测组标签与 pool 完全重合** | **不可归因**；只作逐样本描述，不做 inferential DE |
| 探针级信号（probe-ratio） | probe × cell（`sample_raw_probe_bc_matrix.h5`，54,580 探针） | 多探针基因才可用探针间比值；**MT 基因 1 探针，不可用** | 报告每基因探针数、被排除探针、half-mapped 比例 | 仅对多探针核基因可作探索；MT 基因此法不成立 |
| 293T 潜在状态异质性 | sample × cell；非独立实验重复 | 逐样本评估 QC、doublet、ambient；**ambient 须按 pool 内多样本处理** | 报告每样本保留/排除数与原因 | 可做探索性细胞状态描述；不可泛化为细胞类型 |
| pathway / gene-set 变化 | 可识别 contrast 下的 sample-level pseudobulk ranked genes | 仅基于可解释 contrast；不把公司 cluster-wise DE 当样本级推断 | 报告背景集、multiple-testing、方向一致性及未通过的 contrasts | 仅 pool 内可识别对比可作探索性通路解释 |

后续每项过滤须给出前后 barcode/record 数量、原因及 machine-readable excluded audit。**若某选择改变样本归组、观察单位或对比定义，先确认再运行。**

## 4. 解除阻塞后的工作顺序

1. **先过可行性闸门**：拿到靶点坐标后，查是否落在 `chrM_probe_footprints.csv` 的 600 bp 内（一条命令，10 秒）。结果决定后面有没有间接路线这一支。
2. 补齐 8 样本的真实组别、靶点、control 定义、独立培养/克隆/技术重复关系，以及两 pool 的湿实验关系；更新唯一 SSOT [scripts/0_sample_table.csv](../scripts/0_sample_table.csv)。
3. 确认本轮目标是表达响应，还是另有可读取原始 DNA/RNA 碱基的编辑验证数据。**本批 FRP 不再安排碱基编辑检测分支。**
4. 表达分析从公司 per-sample filtered/raw 矩阵与 QC 起步：逐样本表达 QC → **深度均衡化（必做，因 2.98× 差异）** → ambient 背景（按 pool 内多样本）→ 293T doublet 敏感性。不把常规 STARsolo 当 FRP 重计数器。
5. 构造 **pool 内** sample-level pseudobulk 对比；edited-vs-control 方向仅作逐样本描述。对任何涉及 CBE-Control 的结果做深度敏感性分析。
6. 按研究问题决定是否进行细胞状态与 pathway 分析；出图遵循 V21 palette 与科研作图规范，逐图核对样本点数、标签、可读性。
7. 完成英文研究报告与 `custom_research_report_YYYYMMDD/`，逐条审计主张、来源、过滤影响与结论边界；**交付保持内部审阅、未发送**。

## 5. 继续正式分析所需的信息（请勿仅依据 N1/N2 名称猜测）

| # | 需确认事项 | 为什么这条会改变分析 |
| :--- | :---: | :---: |
| 1 | **8 样本的独立实验单位**：各为独立培养、独立克隆还是技术重复？N1/N2 代表靶点、克隆还是重复？ | 决定 n 与观察单位；决定 pool 内对比是否有重复支撑 |
| 2 | **靶点精确坐标**：若为 mtDNA，给 chrM 1-based 位置 + ref/alt；若为核编码基因，给基因名 | **直接决定靶基因是否被面板定量**（§2.2） |
| 3 | **4 个 edited 的构建信息**：编辑器、gRNA、预期编辑窗口 | 决定间接路线是否成立及窗口范围 |
| 4 | **4 个 controls 的定义与配对**：N1-Control/N2-Control/ABE-Control/CBE-Control 各含什么构建体？哪些样本应互相比？ | 决定可构造哪些 contrast |
| 5 | **两 pool 的实验定义**：是否不同 GEM/建库/处理批次？有无跨 pool 桥接样本？ | 决定 edited-vs-control 是否永久不可识别 |
| 6 | **结论目标**：本轮是否只做 FRP 表达响应？若要判断编辑是否成功，请给已有 DNA amplicon/genotyping 数据的精确路径 | FRP 本身不提供碱基编辑读数（§2.1） |

**关于第 1 条的现有线索（我的推断，待确认，不作为分析依据）：** 客户来信的 **4 edited + 2 base-editor control + 2 gRNA control** 这一计数结构，与 GENEWIZ 的 8 个标签在数量上可完全对应 —— `{N1,N2}×{ABE,CBE}` 为 4 个 edited，`ABE-Control/CBE-Control` 为 2 个 editor-only，`N1-Control/N2-Control` 为 2 个 gRNA-only。**若该读法成立，则 N1/N2 索引的是 gRNA/靶点而非重复，意味着 8 个样本是 8 个互不相同的条件、每条件 n=1，没有任何生物学重复。** 这比"每类对照只有 2 个"更严格，会直接决定能否做任何 inferential 统计。**该读法来自两个独立来源的计数吻合，但仍是推断，必须由客户确认。**

## 6. 【安全提示】客户材料内含明文凭据且已进入 git 历史

实测：[client_materials/1_initial_introduction.txt](../client_materials/1_initial_introduction.txt) 中转录的 GENEWIZ 邮件含 sFTP 主机、用户名与**明文密码**。该文件**已被 git 跟踪并提交**（仓库 `/mnt/ex_8T_SSD/projects_2026Q4`，提交 `aac61f8 p47`），未被 `.gitignore` 覆盖。

**我未对此采取任何处置**（删除、改写历史、改 `.gitignore` 均需用户明确授权；改写历史另受规则禁止）。**需要用户判断的事项：**

| 项目 | 内容 |
| :--- | :---: |
| 精确路径 | `/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/client_materials/1_initial_introduction.txt` |
| 数量 | 1 个文件，1 组凭据（host + user + password），位于文件第 2 段（GENEWIZ 转发邮件部分） |
| 看什么 | 该凭据是否仍有效；GENEWIZ 邮件称数据保留 30 天（2026-10-01 起算，现已过期），凭据可能已失效 |
| 看完做什么 | 决定三件事：① 是否请客户/GENEWIZ 轮换凭据；② 是否把该段凭据从工作区文件中移除（保留邮件其余内容）；③ 是否将 `client_materials/` 加入 `.gitignore`。**git 历史中的那一份只能由你决定如何处理。** 任一项获授权后我再执行 |

## 7. 本文与 Codex 版的差异清单

**Codex 版经复核成立、本文保留的判断：** Flex/FRP 平台识别；碱基编辑不可判定；edited-vs-control 与 pool 完全混杂；不对 Flex 用 STARsolo；本地无 BAM；decision-table 写法；向客户提问的框架。**这些判断是正确的，本文在其基础上增补。**

| # | 类型 | 内容 | 位置 |
| :--- | :---: | :---: | :---: |
| 1 | **新增** | 回答了"8 样↔4 对 FASTQ"这一用户实际疑问，并从原始 reads 实测验证 8 个 probe barcode 全部在位（~96% reads 归属） | §0.1、§1.2 |
| 2 | **证据升级** | 把"读到的是探针"从引用官方文档升级为本地实测（R2 前 50 nt **88.5% 逐字精确命中**，反向互补 0） | §2.1 |
| 3 | **新增（改变可行性判断）** | 探针面板对 chrM 仅覆盖 **3.62%**；mt-tRNA/mt-rRNA 探针 **0** 条；MT 基因 **1 探针/基因**（无基因内对照）；6 个经典致病位点 **0/6** 被覆盖 ⇒ 连间接路线也基本阻断 | §2.2 |
| 4 | **实质修正** | Codex 版把表达差异分析整体判为"不能归因"，**就主方向正确但结论过宽**：pool 内对比（如 ABE vs CBE）不受 pool 混杂影响，是本轮唯一可识别的差异分析入口 | §2.3 |
| 5 | **实质修正** | Codex 版写"每类对照只有 2 个"；按客户 4+2+2 与标签结构的计数吻合，更可能是 **8 个条件各 n=1、零生物学重复**（仍属推断，待确认） | §5 |
| 6 | **新增** | 样本间 read 分配差 **2.98×**；且复杂度最低样本（CBE-Control，26,250 UMI/cell）**正是对照组成员** ⇒ 本项目最可能产生假阳性结论的单点 | §2.4 |
| 7 | **新增** | 同 pool 内 GEM 重叠 **11.4%–17.4%**（`frp_gem_barcode_overlap.csv`）⇒ 同 pool 样本非完全独立，ambient 校正须按多样本处理 | §2.4(c) |
| 8 | **新增资产** | 发现并登记 `sample_raw_probe_bc_matrix.h5`（54,580 探针 × barcode）—— probe-ratio 分析的唯一入口，Codex 版 decision table 未纳入 | §3 |
| 9 | **小更正** | `00_fastq/` 的 8 个 `.md5` 实位于 `md5/` 子目录（Codex 版"16 个文件"数量正确，位置未细分） | §1.1 |
| 10 | **安全提示** | 客户材料含明文 sFTP 凭据且已进入 git 历史，未处置，待用户决定 | §6 |

## 8. 当前完成状态

- **已完成**：读取客户介绍；数据完整性复核（380/380、md5 8/8）；**10x Flex 数据格式实测与存档**（[data_format/](../data_format/) 目录：1 份技术文档 + 1 份机器可读 chrM 足迹表 + 4 个可复现证据脚本）；8 样本 probe barcode 实测验证；探针面板 chrM 覆盖率测定；本计划复核版。
- **尚未完成**：样本设计澄清（§5 的 6 项）、可识别的表达分析、研究报告、交付包。
- **正式下游分析未启动，且这不是阴性分析结果** —— 阻塞项是设计信息，不是技术能力。
- **不对外发送任何材料。** 后续如生成报告或交付包，默认仍处内部审阅阶段。

## 9. 本文事实来源与复核方式

| 本文主张 | 来源 | 复核命令 |
| :--- | :---: | :---: |
| 380/380、md5 8/8 | 下载/看门狗日志 | `tail -1 /mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/logs/download_watch.log` |
| R2 前 50 nt 88.5% 精确命中 | 本地实测 | `conda run -n regular_bioinfo python3 -I /home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/data_format/evidence_scripts/1_r2_read_anatomy.py` |
| 每 pool 4 个 probe BC、~96% 归属 | 本地实测 | `conda run -n regular_bioinfo python3 -I /home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/data_format/evidence_scripts/2_probe_barcode_demux_check.py` |
| chrM 覆盖 3.62%、12 条足迹 | 本地比对 GRCh38 chrM（实测 16,569 bp） | `conda run -n regular_bioinfo python3 -I /home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/data_format/evidence_scripts/3_mt_probe_footprint.py` |
| mt-tRNA/mt-rRNA 探针 0 条 | 面板实测 | `awk -F, 'NR>6{split($3,a,"|"); if(a[2]~/^MT-T|^MT-RNR/) print}' <probe_set.csv> \| wc -l` |
| GEM 重叠 11.4–17.4% | 公司诊断文件 | `cat /srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/genewiz_30-1371703644_1124/30-1371703644/01_analysis/*/multiplexing_analysis/frp_gem_barcode_overlap.csv` |
| read 深度 2.98×、CBE-Control 最低复杂度 | 公司 `metrics_summary.csv` | `grep -E ",Cells,\|Median \|Number of reads in cells" <per_sample_outs>/*/metrics_summary.csv` |
| 凭据文件已进入 git 历史 | 本地 git | `cd /mnt/ex_8T_SSD/projects_2026Q4 && git log --oneline -1 -- 47_Liu_Yuan_scRNA/client_materials/1_initial_introduction.txt` |

**我的判断与推理（非实测，需确认或在报告中标为推断）：**

- **两 pool = 不同 GEM/建库批次**：由"不同 physical library + 不同 FASTQ + 无跨 pool 样本"推断；湿实验流程未知。这是 §2.3 不可识别性结论的前提。
- **8 条件各 n=1**：由客户 4+2+2 计数与 GENEWIZ 标签结构吻合推断（§5），两来源独立但仍为推断。
- **read 分配不均源于上样量/捕获效率差异**：机制推断；湿实验投入量未知。
- **R2 未精确匹配的 11.5% 主要是测序错误**：依据 Q30 RNA read 96.2% 与 `Reads mapped to probe set` 98.4%，未逐条做容错比对确认。
- **"编辑落在连接点影响更大"**：机制推断，本批数据未验证；§2.2③ 的"1 探针/基因无内参"则是实测结构性事实。
