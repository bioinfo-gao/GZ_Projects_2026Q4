# P44 chr11 locus junction-level analysis with read rescue（Analysis 1，S01–S10）— 研究计划与审阅

- **Plan Date（创建）：** 2026-10-24
- **状态：** 2026-10-24 用户授权**全流程**并授权自动安装所需软件；Q1 = P44 独立于 P42 重做；Q2/Q3/Q4 = 采用默认。
- **样本表（SSOT）：** [`scripts/0_sample_table.csv`](../scripts/0_sample_table.csv)
- **客户要求原文：** [`docs/P41_to_41_request.md`](P41_to_41_request.md)（Xin 2026-10-13）；报价往来：[`docs/P42_email_reply_to_Xin_1014_reduced_scope.txt`](P42_email_reply_to_Xin_1014_reduced_scope.txt)

## 更新记录

- 2026-10-24 — 创建。
- 2026-10-24 — 用户批复 §7：Q1 独立重做（P42 只在最后作对照）；Q2 S11–S18 作 artefact background panel；Q3 只整理已有 Manta/DELLY VCF，SNV 只做 chr11；Q4 §4 阈值照表。授权全流程 + 自动安装软件。新建 conda env `p44_assembly_20261024`（SPAdes）。
- 2026-10-25 — 6 样本试跑（`*_test` 目录，独立于正式输出）发现的证据模型补充，均为**数据驱动的 artefact 判定**，§4 阈值未改：
  1. **Adapter read-through split**：insert 中位 161 bp < read 150 bp，大量 read 读入 TruSeq adapter + index + poly-G 尾；poly-G 尾被比对为 supplementary 到 chr2:32,916,3xx（纯 poly-G 区）。⇒ 3′ segment 起点 ≥ adapter 起点的 split 不是基因组 junction（6 样本里删 1,677 条）。**这就是 P41 "chr2:32.916 Mb 吸附区" 的序列层面成因**（P41 已标为 mapping artefact，P44 给出机制，属细化）。
  2. **Low-complexity**：segment / discordant mate / clip 的 Shannon entropy < 1.2 bits 判为 poly-G 类 artefact。
  3. **3′ clip 内含 adapter**：adapter 起点在 clip 前 10 bp 内 ⇒ adapter；更深 ⇒ 只保留 insert 部分（≥ 20 bp）。adapter 之后的比对 segment 一律不作证据。
  4. **Secondary 必须是"竞争性"placement**：AS ≥ primary AS − 10 才计（与 §4 MARGIN 相同）；否则计入 audit `secondary_not_competitive`。
  5. **Event grouping**：同一批 read-end 同时支持多个只差 partner 位置的 junction（例：S04 chr11:114,605,218 → TTAGGG 端粒重复，partner 可放 chrX/chr4/chr2/chr10 末端）⇒ 合并为一个 event，per-sample 计数取并集，列出 alternative placements；有 alternative placement 即属 §4 的 Ambiguous（partner 定位歧义）。
  6. **§4 增补 "Limited" 一级**：≥ 2 fragments、flank 唯一、不确定区间 ≤ 10 bp，但 < 3 fragments 或仅一类证据——原 §4 表在 Well 与 Ambiguous 之间留了空档，按此填补（不放宽 Well）。
  7. SPAdes 对微小局部 read 集合的 k-mer coverage model 报错（exit 64）⇒ 改 `--sc --cov-cutoff off`，且 read 先在 adapter 处截断（≥ 30 bp）。

---

## 1. 范围（用户决定，2026-10-24）

| 项目 | 决定 | 来源 |
| :--- | :---: | :---: |
| 主任务 | 高质量完成 Analysis 1（junction-level analysis with read rescue） | 用户本轮指示 |
| 样本 | **S01–S10 全部 10 个**（比 10-14 邮件报价的 8 个多 S03、S10） | 用户本轮指示 |
| Analysis 2 | 只做**支撑 Analysis 1 判读**的部分（见 §5 Stage 9） | 用户本轮指示 |
| 工时 / 价格 | 不作为约束 | 用户本轮指示 |
| 再比对 | 不做全基因组再比对；只做局部 rescue reads 的再比对 | 客户邮件 |
| 目标位点 | `chr11:114,604,000–114,625,500`（hg38，1-based closed，21,501 bp） | 客户邮件 |

**与 P42 的关系（待用户确认，§7 Q1）：** P42（`/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis/`）是另一条已运行的线，仍在 tmux `p42_*` 中活动。本计划默认 **P44 独立重做**：不读取 P42 的任何结果作为输入，只在最后一步把 P42 的 C1–C4 作为 comparison annotation 对照；P44 不写入 P42 的任何目录。

P42 审计里已经指出的问题，P44 从设计上避免（依据：`42_.../docs/P42_current_analysis_issue_audit_1022.md` §3）：
- 不以 P41 候选或 P41 read-name 封闭输入（问题 A/B）；
- 候选判据不设当前 ~15× 数据无法达到的门槛（P42 的"60× sensitivity context"让任何事件都无法通过，结论只能写成"全部未过 gate"，对客户没有信息量）；
- 交付按客户原话给出**逐 junction 表**（P42 最终只给了 C1–C4 四个候选）。

## 2. 输入与只读勘察

### 2.1 输入（全部只读）

| 输入 | 路径 | 量 |
| :--- | :---: | :---: |
| BAM（bwa-mem2 2.2.1 `-Y`，GATK MarkDuplicates 标记未删） | `/home/gao/projects_2026M7/41_Vedict_18_wgs/analysis/bam/S{01..10}_*.markdup.bam` | 10 个，合计 **481.99 GiB**（`stat`，2026-10-24 23:09；逐样本见样本表 `bam_GiB`） |
| Reference | `/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta`（含 alt/decoy/HLA/chrEBV；与 BAM 同一份） | — |
| Manta per-sample genome-wide | `41_.../analysis/sv/manta_gw/<S>/results/variants/{diploidSV,candidateSV,candidateSmallIndels}.vcf.gz` | 10 |
| Manta joint region | `41_.../analysis/sv/manta_region_joint_20261007/results/` | 1 |
| DELLY discovery / region re-genotype | `41_.../analysis/sv/delly/{discover,region_geno_20261007}/` | 10 + 1 |
| P41 region events（**只作 comparison annotation**） | `41_.../analysis/region/compare/events_all.tsv`、`events_by_sample_long.tsv` | 4,266 events，1,181 in-target |
| P41 depth（100 bp tiles / 10 kb） | `41_.../analysis/coverage/` | — |

### 2.2 邮件中 "150 candidate junctions" 的口径复核

10-14 邮件逐样本数（21/14/18/12/22/24/13/26）与 P41 规则 "in_target 且该样本 `split_nondup_mapq20 ≥ 1`" 基本一致（重算：20/14/17/12/22/24/13/26，合计 148）。扩到 S01–S10：

| 规则 | S01 | S02 | S03 | S04 | S05 | S06 | S07 | S08 | S09 | S10 | 去重后事件数 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `split_nondup_mapq20 ≥ 1` | 20 | 14 | 13 | 17 | 12 | 22 | 24 | 13 | 26 | 18 | **179** |
| `split_nondup ≥ 1`（含 MAPQ<20） | 61 | 77 | 81 | 45 | 70 | 43 | 74 | 46 | 88 | 106 | **672** |

复核脚本：`/tmp/claude-1001/-home-gao--claude/5253fed0-e021-4d8b-b988-67860909f347/scratchpad/cnt.py`（scratchpad，可丢弃）。**这两个数只用于和 P44 的新 discovery 对照，不作为 P44 的输入集合。**

### 2.3 只读勘察（S04，目标区，2026-10-24）

| 记录类别 | n | dup | MAPQ0 | soft-clip | hard-clip | 带 XA | 带 SA |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| primary | 3,561 | 56 | 1 | 1,824 | 0 | 205 | 48 |
| supplementary | 14 | 0 | 7 | 14 | 0 | 8 | 14 |
| unmapped（随 mate 定位） | 4 | 0 | 4 | 0 | 0 | 0 | 0 |
| secondary（0x100） | **0** | — | — | — | — | — | — |

由此得出的设计约束：
1. **BAM 里没有 secondary record**（比对时没用 `-a`）。客户所说的 "multi-mapping / secondary" 证据在本数据里**只以 `XA` tag 形式存在**。⇒ "本属于该位点、却被放到别处"的 reads，只能靠**全基因组扫描 `XA`/`SA`/mate 坐标**找回，区域查询找不到。这是 P44 rescue 的核心增量。
2. **hard-clip = 0**（`-Y`）⇒ 每条 record 的 `SEQ` 都是完整 read，sequence completeness 不需回溯 FASTQ（仍逐样本验证，Stage 1 断言）。
3. 目标区内 primary 几乎无 MAPQ0 ⇒ 目标区本身可唯一定位（与 P41 模拟一致），ambiguity 主要在 **partner 端**（如 P41 的 chr2:32.916 Mb 吸附区）和 114,601.5–114,603.5 kb 低复杂度段。
4. primary 中 51% 带 soft-clip ⇒ 多半是未 trim 的 adapter 端和短插入片段；clip 必须按长度 + 序列（adapter 识别）分层，不能直接当 breakpoint 证据。

## 3. 证据模型 decision table（§2A；正式运行前须经用户审阅）

**Data universe**（先定义，不依赖任何旧候选）：
- **Locus extraction window W** = target ±5 kb = `chr11:114,599,000–114,630,500`（覆盖 114,601.5–603.5 kb 低复杂度段和 S04/S05 的深度边界 114.604–114.609 Mb；±5 kb > 最大 insert size，保证 discordant mate 落入）。
- **Rescue set R(sample)** = 以下任一条件命中的全部 record（含 dup、MAPQ0、QC-fail，全部保留并打标签）：
  (a) record 落在 W；(b) mate 坐标（RNEXT/PNEXT）落在 W；(c) `SA` 任一段落在 W；(d) `XA` 任一 alternative hit 落在 W；(e) 双端均 unmapped 的 pair，经与 W 序列 k-mer/局部比对命中。
  对 R 中每个 QNAME 取回**全部** record（含其在别处的 mate 和 supplementary）。

| 研究问题 | 证据单位 | 保留 / 去重规则 | 判定指标 | 允许的结论 | 不允许的结论 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| 位点有哪些 junction（discovery） | alignment record → read-end → fragment | 所有 record 进 audit 表；按 `(sample,QNAME,read_end)` 合并 alignment context；按 `(sample,QNAME)` 计 fragment；dup 保留但不计数 | split / soft-clip（≥20 bp 且非 adapter）/ discordant 在 ±15 bp 内成簇；P41 事件仅作对照列 | "候选 junction"，附全部支持 record | 用 P41 是否出现来判断真假 |
| rescue 前后支持数 | fragment（主）+ read-end（并列） | **before** = P41 口径（primary/supplementary，非 dup，MAPQ≥20，记录落在窗口内）；**after** = R 中所有 fragment 对 junction 模板再比对后的 alt 支持 | 每 junction × 每样本：before / after fragment 数、after read-end 数、ref-supporting fragment 数 | 支持数随 rescue 的变化及其来源（XA / SA / mate / unmapped） | 把 read-end 数当独立分子数 |
| 断点 bp 结构 | assembled contig + split read 序列 | 局部组装（per sample 与 per cell line pooled 两种）；contig 与 split read 两路独立推断，不一致时都报告 | 两侧坐标、方向、microhomology（query overlap）、inserted sequence（query gap）、repeat-shift 不确定区间 | "bp 分辨率断点 ± 不确定区间" | 在 microhomology / 重复平移区间内给单一坐标而不标区间 |
| 是否因 repeat 而 ambiguous | 每个 junction flank（≥30 bp 两侧）的全部 genomic placement | 保留 XA、contig 的 `-a` 多重比对、BLAST 命中，全部制表 | best vs second-best alignment score margin；flank 的 RepeatMasker 注释；同序列在基因组中的拷贝数 | well supported / ambiguous-repeat / weak / not supported（§4） | 以 primary label 或单一 MAPQ 代替真实位置判断 |
| caller concordance | 每 junction × 每样本 × 每 caller | Manta（gw diploid + candidate，joint region）、DELLY（discovery + re-genotype）全部 record 保留 | 两端都在 ±100 bp 内（与 P41 一致）= concordant；单端近邻 = anchor-only；另报最近 call 距离 | Manta-only / DELLY-only / both / neither | caller 一致 = 独立生物学验证（同一数据） |
| clone-specific（支撑 A1） | fragment | 同一模板、同一计数规则对全部 10 个样本再基因分型 | "absent" 需 ref-spanning fragments ≥10（P41 口径，AF≥0.25 漏检 ≤6%）；否则 "underpowered" | +Dox-only / shared / underpowered | 低深度下的 0 = 不存在 |

**Filter-impact audit（每个过滤都出 before/after）：** W 内 record 总数 → 按类别（primary/supp/unmapped/dup/QC-fail/MAPQ0）；R 中每个来源 (a)–(e) 的 QNAME 数与重叠；soft-clip 被判为 adapter 的数目；每个 cluster 的入选/淘汰原因。全部写 machine-readable TSV。

## 4. Junction 分级判据（预注册；运行后不改，改则新版本并归档）

| 级别 | 条件（全部满足） | 含义 |
| :--- | :--- | :---: |
| **Well supported** | 某样本 after-rescue alt fragment ≥3；≥2 类证据（split / discordant / assembled contig 跨越）；两侧 flank 均可唯一定位（alt 模板 AS − 任何竞争 placement AS ≥ 10，且 flank 不在同源拷贝 ≥2 的序列中）；断点不确定区间 ≤ 10 bp | 可作为 PCR 验证候选 |
| **Ambiguous (repeat)** | alt fragment ≥2，但至少一侧 flank 有 AS 差 <10 的竞争 placement，或断点位于 repeat 平移区间 >10 bp | 有 junction 证据，但位置/伙伴端不能唯一确定 |
| **Weak** | alt fragment = 1（不论 MAPQ） | 单分子，仅列出 |
| **Not supported after rescue** | 再比对后 alt fragment = 0（P41 曾有记录） | 原记录不能被重建 |

阈值 3 / 2 / AS 10 / 10 bp 是本计划提出的默认值，**须经用户确认（§7 Q4）**。另外逐样本报告局部深度和检测功效（AF 0.1/0.25/0.5 下 ≥3 fragment 的 Poisson 概率），使"未检出"能读成"未在该 AF 以上检出"。

## 5. 分步流程（scripts/ 下按序编写 → 自审 → 再运行）

| Stage | 脚本 | 内容 | 主要输出（`analysis/<YYYYMMDD>_<stage>/`） | 预计墙钟 |
| :--- | :---: | :--- | :---: | :---: |
| 0 | `0_common.sh`、`0_sample_table.csv` | 路径、线程预算、样本表 | — | — |
| 1 | `1_window_inventory.py` | W 内全部 record 的 audit 表（FLAG 类别、RNAME/POS/MAPQ/CIGAR/SA/XA/clip/SEQ 长度/dup） + 断言 hard-clip=0、无 secondary | `record_audit.tsv.gz` ×10 | <10 min |
| 2 | `2_genome_scan_rescue.sh` | 每样本流式扫描整个 BAM 一次，按 (b)(c)(d) 收集 QNAME；导出双端 unmapped pairs 供 (e) | `rescue_qnames_by_source.tsv`、`unmapped_pairs.fq.gz` | I/O 限速；482 GiB，10 样本 5 并发约 2–4 h |
| 3 | `3_collect_rescue_reads.sh` | 按 QNAME 取回全部 record（`samtools view -N`）；(e) unmapped pairs 对 W 序列比对筛选；生成每样本 rescue BAM + FASTQ | `rescue.bam`、`rescue_R1/R2.fq.gz` | 1–2 h |
| 4 | `4_evidence_clusters.py` | discovery：split / soft-clip（adapter 剔除）/ discordant 成簇；read-end 与 fragment 两级计数；before 口径计数 | `junction_candidates.tsv`、`filter_impact.tsv` | <1 h |
| 5 | `5_local_assembly.sh` | 每候选 cluster 局部组装（per sample + per cell line pooled）；SvABA（已有，P42 环境）区域模式作为第二组装器；contig `minimap2`/`bwa-mem2 -a` 多重比对 + BLAST | `contigs.fa`、`contig_alignments.tsv` | 1–3 h |
| 6 | `6_junction_genotyping.py` | 为每个 junction 构建 alt / ref 模板；10 个样本的全部 rescue reads 竞争比对；逐样本 after-rescue fragment/read-end 计数、ref 支持数、AF 估计与 95% CI | `junction_by_sample.tsv` | <1 h |
| 7 | `7_breakpoint_annotation.py` | bp 坐标、方向、microhomology、inserted seq、repeat-shift 区间；flank RepeatMasker 注释（UCSC hg38 rmsk 表）与竞争 placement；按 §4 分级 | `junction_table.tsv`（**主交付表**） | <1 h |
| 8 | `8_caller_concordance.py` | 逐 junction × 样本对 Manta / DELLY（及 SvABA）；concordant / anchor-only / none | 合并入 `junction_table.tsv` 的 caller 列 | <30 min |
| 9 | `9_*`（Analysis 2 支撑部分，见下） | | | |
| 10 | `10_manual_curation/` | Well supported + Ambiguous 的逐 junction 人工审阅：IGV（原 BAM + rescue BAM）、contig 比对图、断点序列图；写 per-junction 审阅记录 | `curation/<junction_id>.md` + PNG | 视数量，每个 0.5–1 h |
| 11 | `11_figures.py`、`12_build_delivery.sh` | 图、报告源 `docs/_client_report_source.md`、交付目录 `custom_research_report_YYYYMMDD/` | | |

**Stage 9：Analysis 2 中支撑 Analysis 1 的部分（默认做 9a–9c，9d 待确认）**

| 子项 | 内容 | 为什么支撑 A1 |
| :--- | :--- | :--- |
| 9a 位点 clone-specific 比较 | 用 Stage 6 的统一基因分型，S04/S05 vs S01–S03、S07–S10 vs S06；标 +Dox-only / shared / underpowered，写明 tier、fragment 数、AF 阈值 | 客户问的"只在 +Dox 克隆出现"在该位点上的直接答案 |
| 9b chr11 BAF / LOH（S01–S10） | chr11 联合 germline SNV calling（bcftools mpileup/call 或 GATK HC，限 chr11）；het SNP BAF 沿 chr11q 作图 | 判断 S04/S05 深度降一半是否 = 失去一个 haplotype；断点附近的 het SNP 可把 junction reads 归到丢失 / 保留的 haplotype |
| 9c 深度边界精细化（S04/S05） | 100 bp tile 深度分段 + 置信区间 | 给出 junction 必须落入的区间，作为 discovery 的独立约束（深度边界内没有 junction = 结论本身） |
| 9d（待确认）genome-wide 交付文件 | S01–S10 的 Manta / DELLY unfiltered+filtered VCF 整理、genome-wide SNV/indel VCF（10 样本联合 calling，约 10–20 h 计算）、带 caller 列的 genome-wide per-event TSV | 客户 Deliverables 清单里列了，但属 Analysis 2 的 genome-wide 部分 |

**资源：** 持续 ≤28 核 / 56 线程；Stage 2 是 I/O 瓶颈，限 5 个样本并发、每个 `samtools view -@ 4`。所有长任务在 tmux（`p44_<stage>`）。只写 P44 下新的日期子目录，不覆盖任何已有产物；中途打断可从最后完成的样本续跑（每样本 `.done` 标记）。

**需要安装（运行前在对话中逐项列出再装）：** SPAdes（局部组装；regular_bioinfo 中没有）。RepeatMasker 注释用 UCSC `rmsk.txt.gz` 数据文件，不装软件。SvABA / GRIDSS / BLAST 用已有的 `p42_sv_validation_20261023` 环境（只读调用，不改该环境）。

## 6. 预计交付（客户）

- `junction_table.tsv`：每 junction 一行——bp 断点两端、方向、类型、microhomology、inserted seq、不确定区间、repeat 注释、分级、Manta/DELLY 列、逐样本 before/after fragment 数与 AF。
- `junction_by_sample_long.tsv`、`filter_impact.tsv`、`rescue_source_summary.tsv`（各 rescue 来源贡献多少 reads）。
- 逐 junction 审阅图（Well supported + Ambiguous）、`contigs.fa`、每样本 locus rescue BAM。
- 9a–9c 的表与图；9d 视 §7 Q3 决定。
- 报告：P41 结构 §1–§9，英文；P41 结论不被否定，只在其上补充。

## 7. 开跑前需要用户拍板的问题

| # | 问题 | 默认（推荐） | 影响 |
| :--- | :--- | :---: | :---: |
| Q1 | P44 与 P42 的关系 | **独立重做**，最后用 P42 C1–C4 作对照 | 决定是否复用 P42 中间产物 |
| Q2 | S11–S18 是否作为 **artefact background panel**（只用于判断"多个无关细胞系都出现 = mapping artefact"，不逐样本报告结果） | **是**；Stage 2–3 多扫 8 个 BAM（约 +4 h I/O） | 改变 Ambiguous / artefact 判定的证据 |
| Q3 | 9d genome-wide 交付文件做不做 | 只整理已有的 Manta/DELLY VCF；SNV/indel genome-wide 不做（9b 只做 chr11） | 计算量和交付范围 |
| Q4 | §4 分级阈值（≥3 fragment / ≥2 类证据 / AS 差 10 / 不确定区间 10 bp） | 按表 | 决定"well supported"的数量 |

另：报告覆盖 S03、S10 超出了 10-14 邮件确认的 8 个样本，是否先告知客户由用户决定；本计划不涉及任何对外联系。
