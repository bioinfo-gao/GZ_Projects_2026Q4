# P42 从 P44 junction read-rescue 补充交付中学习的复盘与改进计划

**建立日期：** 2026-10-25  
**文档性质：** 内部方法复盘、可迁移改进清单与待确认计划；不是客户交付物，不改变任何既有分析结果。  
**审阅对象：** `44_wgs_advanced_analysis_from_P41_by_cc/custom_research_report_supplementary_20261025/`。  
**P42 现行交付记录：** `custom_research_report_20261024/`；该目录及其中 `P42_Vedita_18_WGS_1024.md` 维持冻结，不在本复盘中修改。  
**本次状态：** 已完成静态文档、产物、代码结构和视觉审阅；未重新运行任何 P44/P42 分析，未将 P44 数值导入 P42。

## 1. 执行摘要

P44 是一次高质量、范围明确的 chr11 locus supplementary analysis。它最值得学习的不是把原有候选“说得更肯定”，而是将研究问题改写成一个可检验的 evidence model：从完整 locus read universe 出发，保留 multi-mapping / secondary / supplementary / soft-clipped context，以 fragment 而非 alignment record 计数，并让每个 junction 与其最优 genomic placement 竞争。由此，P44 能区分“新 rescue 后出现的真候选”“单 read 弱证据”和“被 unrearranged repeat copy 同样解释的证据”。

P42 的既有设计也有重要优点：它已从 full-window evidence inventory 出发、保留 alignment context、避免将 caller concordance 当成独立 validation，并对 C1–C4 保持克制。两项目的主要差别是最终的 analytical product：P42 以固定 C1–C4 的 comparative adjudication 为核心；P44 以 read rescue 后的全 locus re-genotyping、junction classification、depth boundary 和 B-allele fraction（BAF）的一致性为核心。后者使读者能从“某几个候选是否成立”进一步读到“该 locus 存在何种结构证据、哪些 read 很可能只是 repeat-copy misplacement”。

**本复盘的主结论：** P44 的方法可以用来增强 P42 的证据分层和可审计性，但不应直接替换 P42 的严格 high-confidence clone-/doxycycline-specific 口径。若开展后续工作，应新增一个有新日期的 supplementary delivery；其结论须并列说明 P42 与 P44 的 evidence universe、计数单位和判据差异，不能将两种阈值混成一个“更正”。

## 2. 审阅范围、证据与未做事项

### 2.1 已核查的材料

| 审阅对象 | 核查内容 | 结论 |
| :--- | :---: | :--- |
| P44 supplementary report | §1–§9 的问题定义、判据、方法、结果、限制与交付索引 | 结构完整；关键主张能定位到表、图或局部 BAM。 |
| P44 machine-readable tables | target/window junction tables、per-sample long table、rescue sources、filter impact、caller concordance、depth/BAF tables | 支持 fragment count、classification、baseline state 与过滤影响追溯。 |
| P44 figures | 4 张 locus overview / summary 图及 4 张 junction views 的 PNG 视觉审阅 | 概览图能支持正文；junction view 仍有可读性改进空间。 |
| P44 pipeline scripts | window inventory、rescue、realignment、evidence clustering、genotyping、breakpoint annotation、caller concordance、depth/BAF、table/figure build | 关键 stages 显示了 assertions、audit outputs 或显式 criteria。 |
| P44 integrity manifest | `md5sums_delivery.txt` | 154 个受 manifest 覆盖的文件均通过 MD5；manifest 文件自身未纳入清单，故交付目录共 155 个文件是预期差异。 |
| P42 current report and rationale | P42 1024 report、P42 rationale/process 说明、C1–C4 结果表与交付结构 | 用于确定可迁移项及避免改写历史结论。 |

### 2.2 本复盘不能主张的事项

- 未重跑 P44 pipeline，故不将本次审阅描述为独立重现或对 BAM 原始记录逐条复算。
- 未验证“criteria fixed before results were examined”这一时间顺序；该项来自 P44 报告的陈述，不应被本复盘升格为独立审计事实。
- 未做 PCR、FISH、long-read 或其他独立 assay；P44 与 P42 都是对同一 short-read WGS data 的 algorithmic / mapping re-analysis。
- 未把 P44 的结果写入、替换或“校正”P42 的任何客户交付文件。

## 3. P44 的主要优点

### 3.1 Evidence universe 定义正确，避免 candidate-gated 盲区

P44 的 extraction window 是 target ±5 kb，并纳入：window 内记录、mate 落入 window 的 read pair、SA-tag segment、XA-tag alternative hit，以及经 window-sequence screening 命中的双端 unmapped pair。每个被 rescue 的 read name 会取回其全部 records，而不是只保留旧 candidate 表内的 read name。旧 P41 candidates 仅在后续作为 comparison annotation。

这比只从旧 caller list 或固定候选集开始更适合 repeat locus：primary placement、MAPQ 和旧 candidate label 都不能单独定义真实 locus。P44 的做法保留了“非独立 molecule 但仍有定位价值”的 secondary / supplementary context，同时避免将其重复计数为多条支持。

具体量级也说明其必要性：S01–S10 每个样本 rescue 了 6,328–13,872 个 read pairs，其中 3,920–7,946 对仅因 XA alternative hit 才进入集合。P42 对 C1–C4 所做的 unmapped-mate rescue 得到 5 个 recovered mates、与 C1–C4 均无交集；这个 P42 阴性结果只能解释“unmapped-mate route 对这四个固定候选无增量”，不能外推成“rescue 对整个 locus 无价值”。

### 3.2 Counting unit 与 filter impact 透明

P44 明确定义 read end、fragment 和 event 的不同角色：read end 是 sequence/alignment evidence；fragment（一个 read pair）是 DNA molecule 的计数单位；alternative alignments 保留为 positional evidence，不增加 fragment count。PCR/optical duplicate 仍保留在文件中，但从 support count 排除。

每项重要过滤均可追溯到 `tables/evidence_filter_impact.tsv`，并包含 adapter read-through、poly-G low-complexity tail 与低于 best placement 10 AS points 的 alternative alignment。报告还直接给出 artefact 的序列机制，而非仅写“filtered”。这正是 P42 可进一步采用的范式：每个“没有支持”的结论必须同时保留被排除的 record、排除原因和数量。

### 3.3 用 best-placement competition 取代单一 MAPQ 或 best-hit 判断

P44 的 after-rescue re-genotyping 不以“MAPQ ≥20”作为真定位的充分条件。一个 read end 要跨 junction 两边各至少 20 bp，且 junction-template alignment score 必须比 genome 中最佳 placement 高至少 10 AS points；discordant pair 还须符合两端方向和 fragment-size 约束。随后对 530 个 target-region junction 按 fragment 数、证据类型、flank uniqueness、breakpoint uncertainty 和 repeat ambiguity 分为 Well supported、Limited、Ambiguous (repeat)、Weak 与 Not supported。

该流程给出很有信息量的负结果：530 个 junction 中 272 个在 read 与 genome best placement 竞争后不再获支持，254 个只有一个 fragment，只有 4 个达到至少两个 fragments；没有任何一个达到 Well supported。它把“候选过多”转化为可审计的质量层级，而不是以一个单一过滤阈值简单删除。

### 3.4 Orthogonal-in-same-data evidence 被正确标为“coherent”，而非独立验证

P44 没有把 re-alignment 或 caller overlap 误称为 independent validation，而是把 junction、copy-depth step 和 chr11 BAF/LOH 放在同一个结构解释中：

| P44 观察 | 所支持的有限结论 | 不能推出的结论 |
| :--- | :---: | :--- |
| S04 chr11–telomeric-repeat junction：4 fragments、7 read ends | junction-like evidence 与 S04 depth-step interval 一致 | telomeric partner 的唯一染色体来源，或已由独立实验确认的 terminal deletion |
| S05 fold-back inversion：2 fragments；chr11–chr21 junction：3 fragments | 两个 limited junction 都落在/靠近 S05 depth-step evidence | doxycycline 造成该结构，或单靠 short reads 完整还原 complex architecture |
| distal chr11q BAF/LOH：S04/S05 约 98% baseline-heterozygous sites homozygous，8,004 shared LOH sites 中 99.78% 保留同一 allele | 两个 clone 的 distal loss 与同一 homolog 一致 | 两者为同一 biological clone、同一机制，或精确确定所有 proximal rearrangements |

这种写法是可学习的：同一 WGS 内多类信号的 concordance 可以提高结构解释的连贯性，但仍必须与 independent biological confirmation 区分。

### 3.5 结论强度与证据强度匹配

P44 把 2–4 fragments 的事件定为 `Limited` 或 `Ambiguous (partner end)`，没有因其同时出现 depth/BAF pattern 就提升为 Well supported。报告在 §8.3 又明确：重新比对同一数据不是独立生物学测试；PCR 才是直接 confirmation。这一点尤其重要，因为 P44 的局部 assembly 虽已运行，但四个 ≥2-fragment junction 的 `n_assembled_contigs` 均为 0，因此最终 claims 并未假装得到 contig 级支持。

### 3.6 Delivery 结构利于复核，而不是仅“文件很多”

P44 的 155-file package 按功能分层：报告顶层；`tables/` 供机器读取；`figures/` 与 per-junction views；`locus_BAM/` 供 IGV 检查；`sequences/` 供 primer design；`supporting_analyses/` 供 depth/BAF 复核；`genome_wide_SV_VCF/` 供 caller record 追溯；MD5 manifest 供完整性检验。

对用户而言，尤其有价值的是交付了 **小而可浏览的 locus-specific BAM**，而非要求审阅者先打开数十 GB 的全基因组 BAM 才能查看一个 junction。P42 交付完整 18 个 BAM 的能力仍有必要且不应删除，但可增加类似 P44 的轻量复核层。

### 3.7 客户报告的叙事结构优于技术清单式陈述

P44 的 §2 以 sample group 为单位给出“结果—这意味着什么”，§3 再按主题列 measured findings，§8.2 为每个实际结论给出 section、table、figure 和 BAM 路径。这个层次避免了简单重复：前者服务决策，后者服务复核。

相较之下，P42 的 §2 六条 Core Conclusions 技术上准确，但主要按分析问题排列；P44 的 sample-oriented conclusion table 更容易让客户先理解 S04、S05、controls 与 GM12878 各自意味着什么。

## 4. P44 的缺点、风险与可改进项

### 4.1 科学限制：高质量 short-read re-analysis 仍不是独立 validation

P44 的最重要限制不是流程缺陷，而是 evidence layer 本身：所有 junction、depth 与 BAF 结论来自同一批 WGS。局部 re-alignment、template genotyping 和 caller concordance 可以检验 algorithmic robustness，却不能取代 independent molecule / independent assay。P44 已在 §8.3 正确声明这一点；后续文档也必须继续保留，而不能因证据链更完整就淡化。

### 4.2 有限 fragment support 不足以支撑过度确定性

三个 +Dox-associated junction 各有 2–4 fragments；S04 telomeric partner 本身是多位置等价 sequence。这些结果适合用于 PCR/FISH design 和优先级排序，不适合声称完整结构已经唯一确定。特别是 S04，应始终使用“chr11 end joined to telomeric repeat”而不是“joined to chrX”，因为后者不是数据所能唯一支持的表述。

### 4.3 Figure 5 的标题选择可能造成 partner 唯一性的视觉暗示

`P44_Fig5_P44J0332.png` 顶部以一个 chrX coordinate 表示 B-side，尽管图标题及正文已注明 telomeric repeat 在 chrX、chr2、chr4 和 chr10 interstitial tract 有等价 placement。这可能使快速读图者把“example alignment”误解为“unique chrX partner”。建议未来图中改为：`B side: telomeric-repeat sequence; representative placement shown`，并在图内直接给出 `ambiguous partner` 标签。

### 4.4 Junction view 可增加 read-level 可审计信息

P44 junction views 已展示每个 supporting read end 和 junction sequence，但仍可增加 QNAME anonymized ID、R1/R2、strand、original versus re-aligned record、duplicate state，以及 fragment grouping。这样读者能在不打开 IGV 前先判断“七个 read ends 对应四个 fragments”的去重逻辑；IGV BAM 仍保留为完整上下文的唯一来源。

### 4.5 增量脚本的冻结组织应进一步加强

P44 脚本集中在项目根的 `scripts/`，而非本轮专属的 `scripts/20261025/`。这不影响现有计算，但未来若继续补充分析，容易难以一眼区分“哪一版脚本产出哪一份 delivery”。下一轮应复制所需脚本进入新日期目录，只读旧输入、只写新结果；不得原地修改产出过已交付结果的脚本。

### 4.6 图件总体可读，但概览图可优化信息密度

Figure 1、2、4 的科学信息充分，PDF 缩放后可读；不过 Figure 1/4 的 canvas 很宽，报告内嵌时读者需要较大屏幕或缩放才能同时看清 legend、axis 和细小事件。下轮可保留高分辨率 PDF，同时增加一张更窄的 client-summary panel，或将 detailed panel 定位为 supplementary figure。

## 5. P42 与 P44：哪些可比，哪些不能直接合并

### 5.1 共同基础

两项目均使用 GRCh38-aligned short-read WGS，并面对同一类 chr11 repeat-related placement ambiguity。二者都认可：alignment record、read end、fragment 与 biological event 不能混同；caller agreement 不等于 independent validation；comparator zero 不能自动写成 biological absence。

### 5.2 C2 / EV0623 的关键差异应如实记录

P42 将 P41 EV0623 crosswalk 为 C2，并在既定 standard-eligible rule 下报告 S05 raw support 3、standard support 0、controls raw 0；因此结论为“优先 carry-forward hypothesis，clone/doxycycline specificity 未建立”。P44 的 P44J0289 采用新的 rescue + re-genotyping rule，报告 S05 3 after-rescue fragments、S01–S03 均 0，同时 controls 有 14–18 reference fragments，S04 有 7 reference fragments而被标为 underpowered。P44 对此事件的 classification 是 `Limited`，而不是 Well supported。

这两组结论不应被写成彼此矛盾：它们采用的 evidence universe 和 support criterion 不同。P44 提供了更强的 **case-specific junction candidate** 证据；P42 的 high-confidence clone-/doxycycline-specific gate 还要求更严格的 standard eligibility、local context、sequence-level consistency 和独立 confirmation。因此，P44 不能自动反向宣布 P42 旧结论错误；正确框架是“后续分析在更完整 rescue universe 中将 C2 推进为 limited candidate，仍未建立 doxycycline causal effect 或独立确认”。

| 问题 | P42 1024 的现行口径 | P44 1025 的新增口径 | 后续整合时的正确表述 |
| :--- | :---: | :---: | :--- |
| 数据 universe | full-window inventory 后固定 C1–C4 comparison；rescue 增量聚焦 unmapped mates | 以 window、mate、SA、XA、双 unmapped routes rescue 后重建 locus evidence | 明示 P44 增加了 XA/SA-related mapping evidence；不把 P42 的 5 个 unmapped-mate 结果外推至全部 rescue。 |
| C2/EV0623 support | S05 raw 3、standard 0；case/control zero 不作 absence call | S05 3 fragments；S01–S03 0；S04 underpowered；`Limited` | “S05 的 limited, case-specific candidate evidence strengthened；not an independently confirmed SV.” |
| caller evidence | Manta/DELLY exact two-end concordance 0/20 | target-region 530 events 中 DELLY-only 1、Manta-only 1、both 0 | caller nonconcordance 是 algorithmic context，既不证明也不否定 junction。 |
| biological interpretation | high-confidence clone-/doxycycline-specific rearrangement 未建立 | S04/S05 各有与 depth/BAF pattern 连贯的结构候选 | 使用“structurally coherent but limited short-read evidence”；不写 causal or definitive architecture。 |

## 6. P42 可采用的改进包

### 6.1 优先级 A：扩展 read-rescue 与 re-genotyping evidence model

**目标：** 不再只问“5 个 unmapped mates 是否增加 C1–C4”，而是回答“完整 chr11 window 的 rescued read universe 产生了哪些 junction，哪些在 genome-wide placement competition 后仍被支持”。

**可复用内容：** P44 的 source-route accounting、all-record retrieval、`mem -a` realignment、junction/reference templates、fragment-level count、before/after rescue table、filter-impact audit。

**P42 预期新增输出：**

| 输出 | 最低字段 | 用途 |
| :--- | :---: | :--- |
| `junction_table_target.tsv` | breakends、uncertainty、microhomology、fragment counts、classification、ambiguity、flank placements、per-sample state | 每个 target junction 的主审计表。 |
| `junction_by_sample_long.tsv` | junction × sample、before/after support、reference fragments、state、detection context | 防止将一个 sample 的 0 support 错写成高灵敏度 absence。 |
| `rescue_source_summary.tsv` | 每样本各 route 的 QNAME count、complete/incomplete count、screened unmapped count | 将“rescue 有无价值”限定到具体 route。 |
| `evidence_filter_impact.tsv` | 各过滤原因的 record/read-name/fragment count | 保留负证据和排除理由。 |
| `caller_concordance_by_sample.tsv` | event × sample × callset × match category × record ID | 保留 anchor-only/nonconcordant records。 |

### 6.2 优先级 A：将 depth boundary 与 BAF/LOH 纳入 pre-specified structural coherence layer

**目标：** 在不提高 biological-evidence 等级的前提下，检验 limited junction 是否与 copy-depth step、BAF/LOH、shared retained homolog 等信号一致。

**必须预先写清：**

1. 每类信号的 observation unit、normalisation、过滤和 confidence interval；
2. junction position 与 boundary interval 的关系如何判定；
3. 该 layer 只能支持 `structural coherence`，不能单独使 event 成为 confirmed SV；
4. S04/S05 等非独立 samples 的 biological relationship 不得由 BAF similarity 过度推断。

### 6.3 优先级 B：以 explicit claim ladder 取代二元“有/无 SV”措辞

建议为 P42 supplemental 预注册以下层级，避免把 P44 的 `Limited` 与 P42 的 high-confidence gate 混为同义词：

| 层级 | 最低证据 | 可用措辞 | 禁用推论 |
| :--- | :---: | :---: | :--- |
| Mapping observation | 1 fragment 或仅 alignment/contig clue | `weak observation` / `not interpreted` | true locus、clone-specific、absence in controls |
| Limited junction candidate | ≥2 fragments，明确 template/genome competition，breakpoint uncertainty 已报告 | `limited junction candidate` | independently confirmed, causal, complete rearrangement |
| Repeat-ambiguous junction | support 存在但一端/两端多处等价 placement | `junction-like evidence with ambiguous placement` | specific partner chromosome/locus |
| Structurally coherent limited candidate | Limited junction + pre-specified depth/BAF consistency | `limited short-read evidence coherent with the local structural pattern` | independent biological validation、doxycycline causal effect |
| Confirmed event | 预设 fragment/uniqueness evidence 加独立 assay | `confirmed` | 无 |

### 6.4 优先级 B：交付从“大文件可获得”提升为“关键问题可直接复核”

P42 不应删除现有 genome-wide BAM/VCF delivery；它们承担原始数据访问和用户请求。建议新增、而非替换以下轻量层：

- 每样本 1 个 `locus_rescue.bam` + index，仅含 rescued QNAME 的全部 original/re-aligned records；
- junction/reference template FASTA 和 local-assembly contig FASTA；
- 每个 ≥2-fragment event 的 PNG/PDF junction view；
- 以 client-relative paths 写的 report §9 文件索引；
- manifest checksum 与 report-reference existence check。

### 6.5 优先级 C：报告和图件以决策路径为中心

P42 的技术附录和 14 张图可保留，但正式 report 可采用 P44 的三层叙事：

1. §2：按 S04、S05、baseline、GM12878 等 sample group 写“结果 / 对客户的含义”；
2. §3：按 rescue、junction class、depth/BAF、caller evidence 写 measurable findings；
3. §8.2：每个结论指向 table、figure、BAM 和限制；
4. 图件优先呈现 target repeat architecture、rescue 前后分类变化、关键 junction read-level view、depth/BAF coherence；将细节更高的 technical panel 置入 supporting material。

## 7. 建议的实施顺序与验收标准（仅计划，尚未获准运行）

### 7.1 开始前必须由用户确认的科学选择

这些选择会改变 P42 的 data universe、evidence unit、candidate definition 或 conclusion layer，不能由“自动执行”授权隐含决定。

| 待确认问题 | 选项 | 对结论和工作量的影响 |
| :--- | :---: | :--- |
| 目标范围 | 仅将 P44 输出作为 P42 的解释性 supplementary integration；或以 P42 为根重新运行 P44-style full-window rescue | 后者改变 P42 的 primary analysis universe，工作量和新发现范围显著增加。 |
| 样本范围 | S01–S10 focused cohort；或保留 S11–S18 background panel；或全 18 个样本逐一 re-genotype | 背景样本可帮助识别 recurring artefact，但扩大计算与交付体积。 |
| C2 的正式 claim | 保留 P42 strict high-confidence gate；并新增 P44-style `Limited` tier | 这是推荐方案，可避免“limited candidate”被误读成 confirmed SV。 |
| 交付形式 | 新日期 supplementary package；或仅内部方法评估 | 新结果不可写回 2026-10-24 已冻结 P42 delivery。 |

### 7.2 获得确认后的阶段化执行计划

| 阶段 | 新目录与预期产物 | 验收条件 |
| :--- | :---: | :--- |
| 0. Pre-run design | 新建 dated script directory、decision table、sample manifest、resource plan | 明确 research question、data universe、unit、filter impact 和 permitted conclusion；不修改 P42 historic scripts/results。 |
| 1. Rescue inventory | per-sample record audit、QNAME source ledger、SEQ/hard-clip audit、all-record retrieval | 每一 source route 的输入/输出数可对账；incomplete records 保留并单独标注。 |
| 2. Realignment and discovery | rescued-read re-alignment、evidence clusters、local assembly、candidate audit | 不用旧 C1–C4 关闭 discovery universe；任意筛除均有 machine-readable reason。 |
| 3. Genotyping and annotation | junction/reference templates、fragment support、flank placements、microhomology、uncertainty | alternative placements 不重复计数；每个 event 有 classification 和 per-sample state。 |
| 4. Structural coherence | depth bins/boundary CI、BAF/LOH summary、caller concordance | junction/figure/table 中同一坐标须交叉断言；coherence 不得标为 independent validation。 |
| 5. Delivery | supplementary report、focused BAM/FASTA/tables/figures/checksums | MD5、report paths、table aggregates、figure labels、sample count 全部通过；旧 delivery 不变。 |

### 7.3 必须自动化的自审断言

- 每个 rescue source 的 QNAME、record、complete/incomplete totals 能与输出 BAM/ledger 对账；
- read end、QNAME/fragment、event 三种计数分别输出，禁止相互替代；
- duplicate / low complexity / adapter / alternative-placement filtering 均产生 excluded/audit table；
- 图中显示的 event count、sample count、legend class count 与相应 TSV 一致；
- 断点视图渲染的 mismatch/overlap 与 sequence table 一致；跨步骤坐标偏移由 absolute coordinate subtraction 计算；
- report 中每个 `N–M` range、总数、百分比均从当前 package 的 TSV 重新计算；
- 每个 report 文件路径存在，checksum manifest 覆盖除 manifest 自身外的全部 delivery files；
- prior P42/P41 conclusion 的引用必须按“进一步分析/加强或细化证据”框架写作；若发现实质错误，才另行启动影响审计和更正流程。

## 8. 推荐的后续表述模板

如果用户批准整合 P44-style evidence，建议在未来 supplementary report 中使用以下逻辑，而不是改变历史交付：

> The additional full-window rescue and re-genotyping analysis strengthens the evidence for a limited S05 junction candidate corresponding to P41 EV0623 / P42 C2. The candidate has three rescued fragments and is not observed in the stated DLD1 baseline samples under the new assay definition. Its repeat-context interpretation, limited support, and lack of independent confirmation remain explicit; this result does not by itself establish a doxycycline-causal or independently confirmed rearrangement.

该模板的目的不是预先写入客户报告，而是固定边界：新增分析可以提高 evidence resolution，却不能越过 independent confirmation 的层级。

## 9. 关键来源导航

| 用途 | 路径 |
| :--- | :--- |
| P44 正式 supplementary report | `/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc/custom_research_report_supplementary_20261025/Chr11_Locus_Junction_Read_Rescue_Report_1025.md` |
| P44 target junction audit table | `/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc/custom_research_report_supplementary_20261025/tables/junction_table_target.tsv` |
| P44 filter-impact audit | `/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc/custom_research_report_supplementary_20261025/tables/evidence_filter_impact.tsv` |
| P44 per-sample evidence table | `/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc/custom_research_report_supplementary_20261025/tables/junction_by_sample_long.tsv` |
| P44 sequence/figure review root | `/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc/custom_research_report_supplementary_20261025/figures/` |
| P42 current client report | `/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis/custom_research_report_20261024/P42_Vedita_18_WGS_1024.md` |
| P42 current rationale/process document | `/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis/docs/P42-updated1/P42_analysis_rationale_and_process_20261024_CN.md` |
| P42 C2/caller/rescue source tables | `/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis/custom_research_report_20261024/results/` |

## 10. Change log

- 2026-10-25 — 初建：完成 P44 静态复盘；记录优势、局限、P42 可迁移项、不可直接合并的结论边界及待确认实施计划。未运行分析，未修改交付物。
