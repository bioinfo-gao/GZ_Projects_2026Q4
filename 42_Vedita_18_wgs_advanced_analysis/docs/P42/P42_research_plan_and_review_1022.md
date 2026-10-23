# P42 全队列重复区 SV / CNV 证据重建 — 当前分析指南、研究计划与审阅记录

- **首次建立：** 2026-10-22
- **本次实质修订：** 2026-10-23
- **状态：** 用户已完成旧版审阅并授权正式重跑；本文件是当前唯一有效的分析 guide、研究计划和审阅记录。
- **范围：** 全部 18 个 WGS 样本（S01–S18）。DLD1 的 S04、S05（+Dox）与 S01（parental）、S02/S03（control clones）构成主要因果比较；其余 13 个样本用于队列背景、技术/重复区可比性和细胞系特异性解释，不能被错误地当作 DLD1 的生物学对照。
- **核心位点：** `chr11:114,604,000-114,625,500`（GRCh38/hg38；1-based closed interval）。发现窗口预注册为其两侧各扩展 50 kb：`chr11:114,554,000-114,675,500`。
- **规范样本表（SSOT）：** `scripts/0_sample_table.csv`。
- **历史版本：** 旧的 S01–S10、P41-candidate-gated、仅重比对工作流及实际产物原样复制在 `superseded_for_comparison/P42_candidate_gated_read_rescue/v1_20261022_incomplete_window_evidence/`；它们保留为可比较原始记录，**不是**本轮主线结果。
- **文档真源：** 原独立的 GitHub 方法调研/外部 caller 验证设计已合并入本文件 §4.2.1–§4.2.3；其 2026-10-23 原文只作为历史比较材料保存，不再是并行计划。

## 1. 科学问题与结论边界

已观察到 S04/S05 的覆盖度低于 DLD1 非 Dox 样本，客户的 FISH 亦提示存在染色体片段异常。这是需要严格检验的线索，而不是 SV 已定位或 FISH 已由 WGS 独立证实的结论。本轮逐项回答：

1. S04/S05 是否在 chr11 靶位点或全基因组存在与缺失、重排或易位一致的 copy-number / breakpoint 证据？
2. 若有 junction，两个端点坐标、方向、microhomology 与 insertion sequence 是什么；其定位是否因重复序列而存在同等合理的替代位置？
3. 该事件在 S01/S02/S03 中是否缺失或显著低频，因而可称为 Dox-clone-specific；还是本来就存在于 parental/control？
4. 覆盖度降低是一段连续 CNV、全样本深度/文库差异、可比对性下降，还是多因素组合？

结论随证据等级而定：只有同时通过分子去重、断点两端证据、替代定位审计和对照比较的事件才称高置信。单端 soft clip、MAPQ0 或 caller 标记只能称候选/模糊。FISH 与 WGS 是独立证据类型，不能相互替代。

## 2. 输入、版本与不变性

- 输入 BAM：P41 mark-duplicate BAM/BAI，`/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs/analysis/bam/`；仅只读查询，绝不全基因组重比对、编辑或复制 BAM。
- 输入 SV：P41 Manta 与 DELLY 的 unfiltered/filtered 调用及索引；本轮重新解析记录和原始 read evidence，不把旧汇总表当真值。
- Reference：P41 的 GATK GRCh38 analysis-set；local-assembly contig 和 rescue reads 保持同一 reference，并保留所有替代 placement。
- 输入量：18 个 BAM、约 846 GiB（操作性估计）。新输出仅写入 `analysis/20261023_fullcohort_repeat_resolve/`，日志仅写入 `logs/20261023_fullcohort_repeat_resolve/`。重跑或参数改变创建新阶段/时间戳子目录，绝不覆盖旧结果。

## 3. 已发现的旧版问题与修正

旧版以 P41 的 179 个候选和已记录 QNAME 为入口，并使用 `read.query_name not in wanted or read.is_secondary`。这不是靶窗口的完整证据模型：遗漏了不在 179 个候选中的 read，也在比较替代位置前直接排除了 FLAG `0x100` secondary alignment。对于重复区“真实位置”的问题，这造成确认偏差。

本轮不再用候选表限制 read universe。discovery window 内全部 **primary、secondary（0x100）、supplementary（0x800）** records 都保留并分层；


分析单位明确分两层：

- **alignment-record 层：** 保留 mapping hypotheses，判断重复区定位歧义；record 数不能作为独立支持数。
- **read-end 层：** 以 `(sample, QNAME, read-end)` 归并同一 sequenced end 的 primary/secondary/supplementary 成员，完整保留其定位 context；这是 read-level 支持数，不能替代独立 template 数。
- **fragment/template 层：** 以 `(sample, QNAME)` 归并 paired-end 的 R1/R2；同一 template 在主支持计数中至多计一次。duplicate flag 作为属性保留，在主计数与敏感性分析分别报告。

每一个 filter 输出 filter-impact audit：输入/输出 record、QNAME/read-end、fragment/template、sequence-available 数与原因。P41 179 行仅作为“新全窗口是否重新发现/支持”的注释，绝不是发现闸门。

## 4. 预注册分析设计

### 4.0 Junction evidence decision table（本轮主线的可审查口径）

| 研究问题                           |                                      观察/证据单位                                      |                                                                      保留与去重规则                                                                      |                                     判定指标                                     |                                                              可允许的结论；若删除会漏掉什么                                                              |
| :--------------------------------- | :--------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------: |
| 是否存在断点的原始 short-read 迹象 | BAM alignment record：soft clip、SA/split、supplementary、discordant mate、unmapped mate | discovery window 内 primary/secondary/supplementary 全保留；同一 QNAME+read-end 仅计为一个 molecule，但每个 record 的 FLAG、MAPQ、CIGAR、SA 与坐标均保留 |       record class、clip 方向/长度、SA endpoint、mate endpoint、序列可用性       | 只能说明某类 mapping/junction 迹象是否存在；删掉 secondary 会丢掉替代位置，删 supplementary 会丢掉 split/chimeric 结构，删 MAPQ0 会掩盖 repeat ambiguity |
| 迹象能否支持一个 chr11 junction    |                  molecule：同一 QNAME+read-end 的完整 alignment context                  |                                         molecule 去重只用于独立支持计数；不丢弃其同一分子的 alternative records                                         |  两端坐标与方向、是否双端锚定、独立 molecule 数、同一 molecule 的所有 placement  |                                   可称 read-level candidate 或 ambiguous candidate；不能把多个 records 夸大成多个分子                                   |
| 是否可给出 base-pair breakpoint    |                assembly contig 或跨 junction 的 sequence-bearing molecule                |                               每样本/候选独立 assembly；全 contig placement 均保留；hard/soft clipping 与 SEQ 完整性另审计                               |      contig 双端锚定、方向、microhomology、inserted sequence、竞争 placement      |                                  只有通过全部检查才报告精确 breakpoint；重比对/assembly 属算法复现，不是独立生物学验证                                  |
| 是否为 Dox-clone-specific          |                    统一候选在 5 个 DLD1 样本的 molecule/coverage 证据                    |                                    候选以双端各 50 bp、SV type/方向一致合并；S01/S02/S03 是生物对照，其他细胞系仅背景                                    | S04/S05 与 controls 的 standard/inclusive molecule 数、AF、局部有效深度及检测下限 |                   只能称 Dox-shared、clone-private、candidate 或 negative-at-available-sensitivity；不能将单 clone 直接称 Dox-induced                   |

每一类保留/排除都会输出 filter-impact audit，依次给出 record、QNAME/read-end、molecule、candidate 与 sequence-available 数。此表的唯一原始数据 universe 是 discovery window 及其 unmapped mates；P41 179 candidates 只作后续注释列。

**计数单位（用户于 2026-10-23 确认）：** `QNAME+read-end` 用于正确合并同一 sequenced read-end 的 primary/secondary/supplementary records；paired-end 的 R1、R2 通常来自同一个 library fragment，不能被称为两个 independent molecules。cluster 表同时输出 `(sample,QNAME)` 的 **independent fragment/template count** 与 `(sample,QNAME,read-end)` 的 **independent read-end count**，所有 raw records 均不删除。高置信候选的 `≥3 independent molecules` 主阈值固定以 fragment/template count 评估；read-end count 作为 read-level 辅助证据并列报告。

### 4.1 全窗口 read rescue 与局部重建

1. 每样本枚举 discovery window 全部 alignment records，并从窗口内 mapped record 追取其 unmapped mate。写出 record 清单、molecule ledger 和 FLAG/MAPQ/sequence-availability 审计。
2. 每个分子保留所有 sequence-bearing records；sequence 缺失及 hard/soft clipping 单列。MAPQ0、secondary、supplementary 不因类别丢弃。
3. 以 soft clips、split/SA、discordant pairs、跨断点 contig 和 caller 端点聚类生成**新候选**；先检查每个 read 的全 mapping 集，再评估对 chr11 的支持。
4. 每样本、每候选独立 local assembly（不同样本 raw reads 不混合）。contig 允许多重 placement 地重新比对 GRCh38；只有双端锚定、方向一致且替代 placement 已审计时才报告 base-pair junction。microhomology/insertion 从 contig/read CIGAR 与局部序列比对导出，保存序列来源。
5. 每候选分别给出 standard high-confidence 与 inclusive-rescue 计数：前者为主规则下独立 molecule；后者显示 MAPQ0/secondary/supplementary/mate rescue 贡献。两者绝不混写为“read count”。

**首轮 assembly 候选分层（2026-10-23）：** `5c_triage_assembly_candidates.py` 对 5b 的每一个 bin 写入 `cluster_candidate_audit.tsv`，不删除低支持或 soft-clip-only evidence。为了使首轮 local assembly 有界，只有以下 bin 进入 `assembly_candidate_manifest.tsv`：（i）`SA/split` 且 `independent_fragment_count ≥3`；（ii）跨染色体或相距 >5 kb 的 discordant pair 且 fragment count ≥3；（iii）supplementary evidence 且 fragment count ≥3；或（iv）距可读 caller endpoint ≤1 kb、且 fragment count ≥2 的任一 evidence。未进入首轮 assembly 的 bin 被明确标为 retained context，而非阴性；DLD1 comparison 与后续结果审阅可把它们重新纳入敏感性分析。

**assembly 输入与 sequence 完整性审计（2026-10-23）：** 同一样本内距离 ≤500 bp、其 ±500-bp local windows 重叠的 priority bin 合并成一个 assembly locus，避免把同一重复背景附近的多个 100-bp bin 错当为独立重建。每个 locus 取 window 内的全部 BAM records，并对 priority-evidence QNAME 查回 mate（含 unmapped mate）；所有 primary/secondary/supplementary/unmapped record 写入 `alignment_record_ledger.tsv.gz`，每次无法经 BAM index 查回 mate 的 QNAME、源 FLAG/坐标与异常原因另写 `mate_lookup_failures.tsv`。FASTQ 仅为 assembler 输入而去重：同一 QNAME/read-end 选择最长非空 stored `SEQ`，并保留原 FLAG、CIGAR、MAPQ、坐标、hard-clip 与 sequence-source 审计；非空 `SEQ` 绝不被解释为完整原始 read。

**无法 index-rescue 的 unmapped mate（2026-10-23）：** 若 coordinate-sorted BAM 的 `mate()` 反馈 mate unmapped，不能把该 mate 写成“不存在”：BAM index 本身不能按坐标定位未比对记录。对每个此类 QNAME，顺序扫描其所属样本的完整 BAM 至 EOF，只检索该 QNAME 的 matching records；恢复到的 unmapped `SEQ` 单列为 `reads_recovered_mates.fastq.gz` 并加入同一 local-locus assembly 的 single-end 输入。全扫描仍未恢复者写为 `not_recovered_after_full_bam_scan`，不被解释为 SV 阴性或无原始证据。

**de novo assembly 的资源与解释边界（2026-10-23）：** 每个 sample/local locus 独立用 MEGAHIT 1.2.9 组装（k=21–141、最短 contig 100 bp），最多 4 个 job 并发、每 job 6 threads，持续预算为 24 threads（低于本机 28 physical-core 上限），每 job 内存上限 12 GB。assembly 采用 inclusive rescue sequence input（MAPQ0、secondary、supplementary、duplicate 状态均不作为删除条件），但 contig 的存在不改变独立 fragment 支持数，后者始终回到 QNAME/template 口径。任何失败或 input不足 locus 写独立 status，不能被默认为阴性。

**contig repeat-aware placement（2026-10-23）：** 汇合所有成功 local assembly 的 contig 后，使用同一 GRCh38 / GATK reference 的 `bwa-mem2 mem -a -Y` 一次性重比对，保留 primary、secondary 与 supplementary contig placements，输出每条 placement 的 FLAG、坐标、CIGAR、MAPQ、AS、NM 与每 contig 的 alternative-locus ledger。primary label 不是“真实位置”的充分证据；后续只在双端锚定、方向与 raw read evidence 一致、且替代 placement 被审计的情况下讨论 base-pair breakpoint。

### 4.2 SV caller 与全基因组结构背景

重新读取每样本 Manta/DELLY raw 与 filtered records，按 breakpoint-end、方向和 SV type 标准化，生成 `Manta_only`、`DELLY_only`、`both_concordant`、`both_nonconcordant`、`neither_caller_but_read_assembly_candidate`。caller concordance 提供独立算法一致性，不是高置信 junction 的充分条件；重复区无 caller 记录也不自动否定 read-level evidence。

#### 4.2.1 GitHub 工作流比较带来的方法约束

本轮调研的目的不是以某一个新 caller 取代 Manta/DELLY，而是核对成熟 short-read WGS SV 项目如何区分 discovery、局部组装、multi-mapping 分配和跨样本再分型。

| 项目                                                            |                                                                        其可验证的工作流要点                                                                        |                                                    P42 采用的部分                                                    |                                    不可错误外推的部分                                    |
| :-------------------------------------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------: |
| [Manta](https://github.com/Illumina/manta)                       |                       全基因组 breakend graph 后按候选进行 assembly/评分；discovery scan 排除 secondary、supplementary 与 unmapped record。                       | 既有 Manta 调用作为独立 caller 注释；如历史配置未证明解除 repeat-edge 限制，则新目录补跑靶向 repeat-unmasked Manta。 |    Manta 未使用的 record 不是无生物学价值；其内部过滤不得变成 P42 原始证据的删除规则。    |
| [GRIDSS2](https://github.com/PapenfussLab/gridss)                | 靶向提取 ROI-overlap fragment 的 mate 与 supplementary record；可报告 single breakend，且将 multi-mapping unique allocation 前的 raw breakpoint 与分配后支持区分。 |           建立独立的 raw-versus-allocated caller evidence 表，保留 single breakend 与 assembly sequence。           | GRIDSS 要求每个 read 一条主要 alignment，不能直接吞入 secondary；这只是 caller 输入约束。 |
| [SvABA](https://github.com/walaj/svaba)                          |                                    对目标 BED/区域 local assembly、BWA-MEM 重比对 contig，输出 VCF 与每样本 detailed evidence。                                    |                                 作为 MEGAHIT + all-placement contig 的独立算法复核。                                 |            caller score 或一条 primary contig placement 不能证明真实唯一位置。            |
| [UCLA call-gSV](https://github.com/uclahs-cds/pipeline-call-gSV) |                                                       将 SV/CNV discovery 与已发现候选的 regenotyping 分开。                                                       |                        DLD1 的 clone-specific 判定在固定候选坐标对 S01–S05 用同一规则重测。                        |                控制样本 VCF 没有一行 call 不能替代可量化的阴性/检测下限。                |
| [nf-core/sarek](https://github.com/nf-core/sarek)                |                                                 以可复现输入、caller 标准化、QC 和可审计交付组织多 caller 工作流。                                                 |                            记录版本、参数、manifest、per-sample QC 和标准化 caller 字段。                            |                   多 caller 的严格交集不能作为重复区候选的硬淘汰门槛。                   |

Manta 的 record 过滤行为以其[方法文档](https://github.com/Illumina/manta/blob/master/docs/methods/primary/methods.tex)为准；GRIDSS 的 targeted extraction、single-breakend、raw-before-allocation 输出与 joint calling 建议以其[README](https://github.com/PapenfussLab/gridss)和[QuickStart](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md)为准；SvABA 的 targeted local assembly 和 per-sample evidence 以其[官方 README](https://github.com/walaj/svaba)为准。

#### 4.2.2 原始证据轨与 caller-compatible 验证轨

P42 明确维持两条不可互相替代的轨道：

1. **full-evidence track：** 18 个样本 discovery window、target-overlap mate 和所有 primary/secondary/supplementary/MAPQ0/unmapped-mate context 均保留在 ledger。它回答“有哪些 mapping hypotheses、哪些位置存在歧义”。
2. **caller-compatible track：** 为每个外部 caller 按其合法输入规则新建衍生输入和 filter-impact audit。它回答“该算法在明确的输入限制下是否也能重现候选”。任何 caller-compatible 排除均写明 QNAME、FLAG、坐标和原因，原始 ledger 不删除。

二者的计数与结论边界如下：

| 研究问题                                    |                             evidence universe / 单位                             |                                                 计数或过滤口径                                                 |                                               允许结论                                               |
| :------------------------------------------ | :-------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------: |
| repeat 区有哪些异常 mapping/junction 线索？ |          18 样本的 record、QNAME+read-end、QNAME/template 三层 ledger。          |                        先完整保留 FLAG 类别；template count 为主，read-end count 并列。                        |                     观察到某类证据、或在可用数据中未观察到；不单独宣布唯一断点。                     |
| 局部 assembly 是否支持候选？                | full-evidence track 的 59 local locus；FASTQ 仅是组装输入，完整 BAM ledger 另存。 |               同 QNAME/read-end FASTQ 选最长可用 stored `SEQ`；其他 alignment context 不删除。               | `assembly-supported`、`unresolved` 或 `repeat-ambiguous`；primary placement 不是充分定位证据。 |
| 外部 caller 是否算法上重现？                |    靶向 caller-compatible 衍生输入；S01–S05 为 DLD1 主比较，18 样本保留背景。    | GRIDSS secondary exclusion 只发生在衍生输入；raw/allocated support、single breakend、assembly score 分开报告。 |      `algorithmically concordant`、`discordant` 或 `not evaluable`；并非独立生物学验证。      |
| 是否 Dox clone-specific？                   |                    固定候选在 S01–S05 的同规则逐样本再分型。                    |    不用 VCF 缺失代替 0；报告 template/read-end support、reference support、coverage、ambiguity 与检测下限。    | `Dox-enriched/clone-specific candidate`、`shared`、`insufficient evidence` 或 `ambiguous`。 |

#### 4.2.3 外部验证的实现与审计状态

GRIDSS 2.13.2 与 SvABA 1.2.0 安装在独立环境 `p42_sv_validation_20261023`，不改变 `regular_bioinfo`。在运行 GRIDSS 前，已对 S04 BAM 与 P41 同一 `Homo_sapiens_assembly38.fasta` 的全 contig name/length dictionary 做逐行比对，结果一致；传统 BWA index 只在 `analysis/20261023_fullcohort_repeat_resolve/06_external_sv_validation/gridss_reference_v1/` 建立，不修改共享 reference 目录。待 index 完成后，依次生成 targeted caller-compatible input、其 filter-impact audit、GRIDSS/SvABA 结果和跨样本 regenotyping 矩阵；每个工具都写入独立目录，绝不覆盖 Manta/DELLY 的历史全基因组 VCF。

已完成的 local-assembly 输入与输出核验：5 个 `mate is unmapped` QNAME 均已通过所属完整 BAM 的顺序扫描回收一个 150-nt R2；59/59 local locus 完成 MEGAHIT、无 failed/insufficient-sequence status，生成 90 条 contig（121,426 bp），并已作 hg38 all-placement alignment。此时尚未从这些输出发布任何 SV 生物学结论，必须先完成 multi-placement/ambiguity 分类。

### 4.3 DLD1 clone-specific 规则

主比较为 S04/S05 对 S01/S02/S03；S06–S18 仅作细胞系背景。候选合并条件为同 SV type、同方向、两个端点各在 50 bp 内；base-pair assembly event 以精确 sequence 优先。每候选列出五个 DLD1 样本的独立 molecule、总可评估 molecule/局部 depth、standard/inclusive AF、record-class 构成及 caller support。

- **Dox-shared high confidence：** S04 和 S05 各 ≥3 independent molecules，至少一个样本有双端/assembly junction；standard AF ≥5%；S01/S02/S03 同 junction 各为 0 或低于预计算检测下限。
- **Dox-private high confidence：** 仅 S04 或 S05 满足上述 molecule/AF/双端条件，全部 DLD1 对照不满足；只称 clone-private，不泛称 Dox-induced。
- **candidate / ambiguous：** 2 molecules、单端证据、MAPQ0/secondary 主导或多个等价 placement；保留、图示、计数，但不宣称定位明确或 clone-specific。
- **negative at available sensitivity：** 对没有观察的样本，给有效覆盖度和检测下限；不写“事件不存在”。

`≥3 molecule / ≥5%` 是初始阈值，原始计数始终并列报告。若 5% 小于 3 molecule，则 3 molecule 仍为下限；局部有效深度不足时降级为不可判定。

### 4.4 Coverage / CNV 验证轨

全队列做两种分辨率 coverage/CNV 审计：全基因组 100-kb bins 检查样本深度和大片段 CN；chr11 10-kb、靶位点 1-kb bins 定位连续 CN boundary。standard CN 使用 primary、non-duplicate、非 QC-fail、非 supplementary、MAPQ≥20 records；另给 MAPQ≥0 inclusive track。每轨明确 flag、MAPQ、normalisation 与 mappability 限制，禁止把重复区 mapping-depth 当真实 CN。

### 4.5 SNV/indel 与交付

对 DLD1 五样本核验同一 reference、同一过滤定义下的 SNV/indel VCF；若已有 P41 callset 的 provenance/完整性不合格，才在新目录补充 calling。全 18 样本 BAM/BAI 和 Manta/DELLY VCF 是既有交付源；本轮生成 manifest、校验和与 annotated per-event TSV，不复制 846-GiB BAM 到 P42 工作区。

## 5. 图片交付（PNG、PDF、来源 TSV 与脚本齐全）

至少制作下列 12 组图，且所有图片均可从保存的 TSV/脚本重建：

1. 全 18 样本 genome-wide normalized coverage heatmap/overview；
2. S04/S05 相对 S01–S03 的 chromosome-level coverage 差异；
3. chr11 10-kb CN tracks（standard / inclusive 并列）；
4. 靶窗口 1-kb depth、MAPQ、有效 molecule tracks；
5. 全样本 FLAG / MAPQ / sequence-availability audit；
6. read rescue 前后 independent-molecule 计数与类别构成；
7. Manta–DELLY–assembly/read evidence caller-concordance 图；
8. DLD1 五样本候选 junction AF/support heatmap；
9. 每个重点 junction 的端点、方向、microhomology 示意；
10. 每重点事件 read/contig alignment evidence panels（含替代 placement/歧义注记）；
11. local-assembly contig 双端锚定与替代位置比较；
12. FISH 观察、coverage、caller、junction evidence 的支持/未支持边界总结图。
13. GRIDSS raw-versus-unique-allocated support、single breakend、SvABA evidence 与 Manta/DELLY 的并列一致性图。

图用色盲友好、白底可读的固定样本配色；不使用 pie chart。若没有 high-confidence junction，图 9–11 仍展示最强候选及不能定位的原因。

## 6. 执行步骤与质量闸门

1. `scripts/20261023/0_common.sh`：路径、reference、区域、全样本断言。
2. `1_inventory_full_window.py`：完整 record/molecule/mate rescue inventory。闸门：18 个完成标志；secondary/supplementary/MAPQ0/SEQ=* 都有显式计数。
3. `2_make_coverage_bins.py` 与 `3_coverage_cnv.sh`：两套 coverage。闸门：每 bin 可追到样本、flag/MAPQ 口径，样本总深度与 BAM 规模相符。
4. `4_parse_callers.py`：Manta/DELLY raw/filtered 标准化和匹配。闸门：明确每份 caller 文件是否存在、可读及 record 数；空文件绝不伪装阴性。
5. `5_cluster_junction_evidence.py`：从全窗口、非 P41 candidate list 聚类断点，并同时写 fragment/template 主支持和 read-end 辅助支持。
6. `5c_triage_assembly_candidates.py`：对全部 bin 写出含选择理由的 audit，并生成首轮 bounded local-assembly manifest；未选中 bin 不删除。
7. `6a_prepare_local_assembly_inputs.py`：将 priority bins 合并为同样本 local loci，写 FASTQ、全 record ledger、mate-rescue 与 sequence-completeness audit。
8. `6b_run_local_assembly.sh`：按样本/local locus 运行 MEGAHIT，并逐 locus 写成功、输入不足或失败状态。
9. `6c_align_local_contigs.sh` 与 `6d_summarize_contig_placements.py`：以 `-a -Y` 全 placement 比对 contig 并输出定位歧义 ledger。
10. `6e_prepare_gridss_reference.sh`：对 BAM/reference dictionary 进行硬性一致性检查，并只在项目私有目录建立 GRIDSS 所需传统 BWA index。
11. `6f+` external validation scripts：建立 targeted caller-compatible inputs 与排除审计；运行 GRIDSS/SvABA，必要时补跑 repeat-unmasked targeted Manta；输出 raw/allocated caller evidence 与 fixed-candidate regenotyping 矩阵。
12. `7_dld1_comparison.py`：DLD1 阈值、AF、检测下限与敏感性分析。
13. `8_make_figures.py` / `9_write_report.py`：图片、报告与主张→来源对照。

超过数分钟的步骤仅在 tmux 运行，单一长任务最多 28 physical cores / 56 threads，并启动 watchdog。先做语法检查、小范围 dry-run、BAM/VCF 可读性检查；失败只新建阶段目录，不清除任何旧文件。

## 7. 风险与停止条件

- 该重复区可能不能用短读长唯一定位；这是有价值的限制结论，不能以挑选 primary placement 掩盖。
- FISH 的片段观察不能直接换算为 WGS breakpoint；最终并列呈现而非合并计数。
- assembly 输入不足、BAM sequence 缺失比例高或 caller 源文件缺失时，量化缺失及对结论影响，继续可执行分支。
- 全窗口和旧 P41 candidate-gated 的事件 universe 若实质冲突，旧版不删除；新旧并列并解释差异来源。

## 8. 变更记录

- 2026-10-23 — 将开源 WGS SV workflow 调研、外部 caller 验证 decision table 和已完成的 assembly/input 完整性审计合并进本文件，使其成为唯一当前分析 guide、plan 和 review；原独立 benchmark 文档归入 `docs/P42/superseded_for_comparison/`，仅供历史逐项比较。
- 2026-10-23 — 用户审阅完成并要求重新制定严谨计划、分析全部 18 样本；旧版 179 candidates/QNAME gate 及直接排除 secondary 的主线废止。旧结果/脚本复制归档；当前改为 full-window、molecule-aware、repeat-aware，新增 CNV、local assembly、DLD1 clone-specific 与大量图形交付。
- 2026-10-23 — 用户确认主支持单位为 `(sample,QNAME)` fragment/template，`(sample,QNAME,read-end)` 仅作并列 read-level 指标；聚类脚本据此更正并恢复运行。
- 2026-10-23 — 5b 输出显示大量 repeat-associated soft-clip bins；新增 5c 作为全量 bin 的可审计分层，而非把单类 soft clip 直接解释为 SV。首轮 assembly 仅运行于高信息量或 caller-anchored evidence，所有未入选 bin 仍保留在审计表。
- 2026-10-23 — 5c 首次运行发现 supplementary cluster-key parser 将两段 key 错误按三段读取，使 16 个 supplementary priority bin 缺失 anchor 坐标；首次脚本与输出归档于 `superseded_for_comparison/P42_assembly_triage/v1_20261023_supplementary_key_parser/`，修复后写入新的 `assembly_triage_v2/`，不覆盖旧产物。
- 2026-10-23 — 在 corrected 5c v2 的 priority manifest 基础上，新增 6a：同样本相邻 priority bin 以重叠的 ±500-bp window 合并，先写全 alignment ledger 与 sequence-completeness audit，再将最长 stored `SEQ` 作为 de novo assembly 输入；FASTQ 去重不删除原始 alignment context。
- 2026-10-23 — 新增 6b：以每 sample/local locus 独立 MEGAHIT assembly，4×6 threads 的资源上限；每个 locus 记录成功、输入不足或失败状态，assembly contig 仅为 mapping/sequence reconstruction evidence，不重定义 fragment-level support。
- 2026-10-23 — 新增 6c/6d：成功 contig 统一以 `bwa-mem2 mem -a -Y` 重比对同一 GRCh38 reference；输出所有 placement 及 alternative-locus ledger，禁止将 primary label 单独解释为真实 locus。
- 2026-10-23 — 6a 全队列初次运行记录 5 次 mate lookup 失败；在启动 assembly 前补写 QNAME/FLAG/坐标/异常原因的 machine-readable failure audit。该 failure 表示无法经 BAM index 取回 mate，不等同于无原始 read 证据或 SV 阴性。
- 2026-10-23 — 5 条 failure 均明确为 `mate is unmapped`。因原始 FASTQ 在 P41 文档所记 Dropbox 路径已不在磁盘，以完整 BAM 顺序扫描作为可用的 source-record rescue：仅扫描受影响的 5 个样本，恢复到的 unmapped mate 追加进入对应 local assembly，逐 QNAME 报告恢复成败。
- 2026-10-23 — mate rescue、local assembly 与 contig placement 以可见的 6ac stage chain 串接：只有 rescue completion marker 存在才可进入 assembly，任何前置阶段失败均终止 chain 并保留已有日志/产物，禁止静默跳步。
- 2026-10-22 — 旧版建立；原文与实际脚本在顶部指定归档中保留，供逐项比较。
