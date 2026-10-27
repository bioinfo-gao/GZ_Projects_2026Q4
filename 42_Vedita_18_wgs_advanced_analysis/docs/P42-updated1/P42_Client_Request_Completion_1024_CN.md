# P42 补充完成报告：针对 P41 第 7.3 节的客户问题（中文版）

**报告日期：** 2026-10-24  
**编写人：** Zhen Gao, PhD，Principal Bioinformatics Scientist，Athenomics  
**分析平台：** Linux HPC server  
**物种与参考：** *Homo sapiens*，GRCh38  
**样本背景：** 依客户请求中命名的 DLD1 parental population 和 DLD1 clone 样本，含/不含 doxycycline；客户未说明组织来源。  
**对应英文原件：** `custom_research_report_supplementary_20261024_v8/P42_Client_Request_Completion_1024.md`  
**状态：** 仅供内部中文阅读与科学审阅；英文原件和本中文版均未发送给客户。

## 1. 目标

本补充分析回答 P41 chr11 repeat-region whole-genome sequencing（WGS）报告第 7.3 节之后提出的问题。分析只使用既有 GRCh38-aligned records；未重新比对，也未新增测序。本轮新增工作包括：

1. 在 read-end 层级解析现有 P42 C2 split-read alignment representation；
2. 记录 Manta 与 DELLY 对该候选是否给出 exact two-end breakend（BND）concordance；
3. 定量已完成的 unmapped-mate rescue 对 C1-C4 的影响；
4. 按预先规定的 DLD1 parental/control 对 +doxycycline 的判定标准进行比较；
5. 盘点客户要求的 genome-wide source files。

英文补充包当时只提供 source-file inventory，并未复制约 0.91 TB 的 BAM；正在构建的单一客户包将提供实际 BAM/index 和 caller 文件。

## 2. 核心结论

1. **C2 的技术描述较 P41 更精确。** P41 的 EV0623 可与 P42 C2 对照：P42 表示为 chr11:114,607,191 与 chr21:20,648,242。3 个不同的 SAM read-name（QNAME）和 5 个已解析 split-read ends 均显示 2-base-pair（bp）的 `AT` query overlap 与 0-bp unaligned query gap。这比 P41 的 split-read 计数提供了更多细节，但仍是 repeat-context alignment evidence，不能视作独立 biological confirmation。
2. **在已审计的 DLD1 records 中，C1-C4 均未获得 Manta 或 DELLY 的 exact two-end concordance。** 20 个 candidate-by-sample comparisons 中，Manta exact match 为 0/20，DELLY 为 0/20，两者同时 exact 为 0/20。C4 有 5 行重复出现的 DELLY chr11 insertion anchor，但没有 chr3 partner，因此被明确标为 anchor-only/nonconcordant，而不是 C4 support。
3. **C1-C4 均未达到 P42 预设的 high-confidence clone-specific 或 doxycycline-specific gate。** C2 仅在 S05 有 3 个 raw dual-endpoint QNAME，standard-eligible fragments 为 0，且 local depth 低于预设的 60x sensitivity context。raw-support/depth ratio 只是描述性 proxy，不是经验证的 allele fraction。
4. **Read rescue 没有增加 C1-C4 的 support。** 全 locus 保留了 5 个 recovered unmapped-mate QNAME，但没有一个与 C1-C4 的 dual-endpoint QNAME 重叠；各候选 rescue 前后 raw count 相同。
5. **客户要求的 source-data inventory 可用于 transfer planning，但有一项明确缺口。** 18 个 BAM 与 18 个 BAM index 均存在（907.44 GB；845.12 GiB），Manta/DELLY source files 已盘点；未发现具 provenance 的 project-wide SNV/small-indel VCF，因此不得将其表述为可用。

第 3 节按主题组织实测发现。第 8.2 节给出每项结论的 evidence path 与解释边界；第 8.3 节说明本数据不能建立的结论；第 8.4 节将下一步科学证据与行政行动分开。

## 3. 关键发现

1. **P41→P42 精确化：** P42 将 P41 EV0623 对照到固定的 C2 geometry；两端各 1 bp 的差异是同一候选的 representation difference，不能作为第二个事件。
2. **Junction structure：** C2 的 split alignment representation 中，read ends 均具有一致的 `AT` overlap，两个 mapped blocks 之间未观察到 unaligned query bases。
3. **Caller comparison：** Manta/DELLY 采用 two-end、50-bp tolerance rule；邻近的 single-anchor insertion 不会被静默提升为 trans-chromosomal BND match。
4. **Clone comparison：** 客户指定的 5 个 DLD1 样本以同一 candidate set、support rule、depth context 与 sequence-level requirement 比较；local depth 偏低，因此 negative count 不能解释为达到既定灵敏度的排除。
5. **交付准备度：** manifest 区分“已有的 source files”与“尚未证实存在的 files”，记录 size/index status，而不是将约 0.91 TB BAM 再复制到当时的补充包。

## 4. 样本信息

| Sample | 客户请求中的 DLD1 角色 | 客户请求中的 treatment label | 本补充中的作用 |
| :--- | :---: | :---: | :--- |
| S01_DLD1_pop | Parental population | 未说明为 +doxycycline | Comparison reference |
| S02_DLD1_c1_ctrl | Control clone | control | Comparison reference |
| S03_DLD1_c2_ctrl | Control clone | control | Comparison reference |
| S04_DLD1_c3_dox | DLD1 clone | +doxycycline | Case group |
| S05_DLD1_c4_dox | DLD1 clone | +doxycycline | Case group |

## 5. 分析 rationale 与判定标准

| 问题 | Observation / evidence unit | 保留与计数规则 | 判定标准 | 允许的解释 |
| :--- | :---: | :---: | :--- |
| C2 在 read-end 层级是否内部一致？ | 按 `(QNAME, read end)` 分组的、带 sequence 的 primary/supplementary alignment records | 保留原始 FLAG、CIGAR、MAPQ、coordinates、query intervals；QNAME 是 observation count，不是 biological replicate | 比较 chr11/chr21 query overlap 或 gap 及 overlap sequence | 与 microhomology 一致的 split-read alignment representation；不等于唯一确定的 breakpoint |
| Manta 与 DELLY 是否支持 C1-C4？ | 72 个 source inputs（18 samples × 4 caller stages）的全部 audited records | 保留所有 source records；缺失 records 不转写为 biological absence | 任意 partner order 下，两个 BND ends 都须在 fixed candidate 的 50 bp 内 | 同一 WGS 数据中的 algorithmic concordance，不是独立 validation |
| DLD1 +doxycycline candidate 是否 high-confidence？ | S01-S05 的 fixed-candidate dual-endpoint counts 与 local depth | raw 与 standard-eligible counts 分开；repeat-context support/depth 不作 allele-fraction estimator | 至少 3 个 standard-eligible independent QNAME fragments、至少 60x local depth context、且有 sequence-level junction | 仅当所有 gates 均通过才可称 high-confidence；不设 direct AF cutoff |
| Unmapped-mate rescue 是否改变 C1-C4？ | Recovered unmapped-mate QNAME 与每个 candidate raw dual-endpoint QNAME 的 intersection | 保留 5 个 recovered mates，且报告包括 zero 的所有 intersections | 每个 candidate 的 rescue 前/后 raw count 与 overlap size | 对 C1-C4 的实测 rescue impact；zero overlap 不代表 rescue 在其他位置无价值 |
| 哪些 requested data types 可用？ | source file、index、existence flag、byte size | inventory 为 read-only；当时的大型 BAM/VCF/BCF 不复制到 scientific supplement | 可用性要求有存在且具 provenance 的 source path；Manta candidate-stage file 不改称“unfiltered” | transfer-planning manifest；必须明确 unavailable categories |

## 6. 方法

分析从已完成的 P42 whole-window evidence inventory、local assembly 与 read-rescue outputs 开始；在存在时保留 MAPQ 0、multi-mapping、secondary、supplementary、soft-clipped 与 unmapped-mate context。本补充不重新比对现有 BAM。

对于 C2，读取既有 raw-support audit 中所有带 sequence 的 alignment records，并在单个 read-end 层级解析。由 CIGAR 推导 query coordinates，再比较 chr11 与 chr21 blocks 的 overlap/gap。这一做法保留 primary 与 supplementary records 作为 structural context，同时避免将多个 alignment records 误称为独立 molecules。

对于 caller concordance，72 个 Manta/DELLY source files 与各 C1-C4 geometry 比较。Exact match 要求两个 BND endpoints 在任意 partner order 下均落于 50 bp tolerance；仅有一个邻近 endpoint 或没有 fixed candidate partner 的 record 仍保留并标记为 nonconcordant。Clone comparison 从已完成 P42 fixed-candidate matrix 读取 raw dual-endpoint counts、standard-eligible counts 与 standard local depth。Data manifest 只记录 availability/size，不改变、复制或发送 source dataset。

## 7. 结果

### 7.1 P41 第 7.3 节到 P42 C2 的精确化

P41 EV0623 报告为 chr11:114,607,192 到 chr21:20,648,241，S05 有 3 split reads、0 discordant pairs。P42 C2 表示为 chr11:114,607,191 到 chr21:20,648,242。两端各 1 bp 的差异是在对照同一 candidate geometry 时的 breakpoint-representation difference。P42 先前已建立 C2 in-silico assay design；本节新增的是对既有 alignment evidence 的独立 read-end accounting。

### 7.2 C2 read-end junction 细节

| C2 测量项 | Observed value | Interpretation boundary |
| :--- | :---: | :--- |
| Distinct QNAMEs | 3 | 作为保留的 observations 计数，不是 independent biological validation |
| Resolved split-read ends | 5 | 全部来自 S05，并保留原始 alignment context |
| chr11/chr21 MAPQ | 每个 resolved block 均为 60 | 高 alignment confidence 本身不能消除 repeat-context placement ambiguity |
| Query overlap | 5/5 read ends 均为 2 bp | 每个 overlap sequence 均为 `AT` |
| Unaligned query gap | 5/5 read ends 均为 0 bp | 未观察到两个 mapped blocks 间有未比对 query sequence |
| 允许的结论 | Alignment representation 与 2-bp microhomology 一致 | 不是经独立验证的、唯一的 biological breakpoint |

### 7.3 Manta 与 DELLY concordance

| Candidate | S01-S05 中 exact Manta two-end matches | S01-S05 中 exact DELLY two-end matches | 两 caller 均 exact | Anchor-only/nonconcordant rows | 解释 |
| :--- | :---: | :---: | :---: | :---: | :--- |
| C1 | 0 | 0 | 0 | 0 | 已审计 source 中无 exact caller match |
| C2 | 0 | 0 | 0 | 0 | 已审计 source 中无 exact caller match |
| C3 | 0 | 0 | 0 | 0 | 已审计 source 中无 exact caller match |
| C4 | 0 | 0 | 0 | 5 | 重复 DELLY chr11 insertion anchor 没有 chr3 partner，不是 C4 concordance |
| Total | 0/20 | 0/20 | 0/20 | 5/20 | 在 two-end rule 下，Manta、DELLY 或两者均未 exact call C1-C4 geometry |

### 7.4 Read-rescue impact

| Candidate | Rescue 前 raw dual-endpoint QNAMEs | Rescued-mate QNAME overlap | Rescue 后 raw dual-endpoint QNAMEs | 结果 |
| :--- | :---: | :---: | :---: | :--- |
| C1 | 1 | 0 | 1 | Unchanged |
| C2 | 3 | 0 | 3 | Unchanged |
| C3 | 2 | 0 | 2 | Unchanged |
| C4 | 1 | 0 | 1 | Unchanged |

全 locus 共 recovered 5 个 unmapped-mate QNAME；没有一个与以上 4 个 candidates 的 dual-endpoint QNAME sets 相交。

### 7.5 DLD1 parental/control 与 +doxycycline 的比较

| Candidate | S04 / S05 raw support | S04 / S05 standard-eligible support | S01-S03 raw support | 适用结论 |
| :--- | :---: | :---: | :---: | :--- |
| C1 | 0 / 1 | 0 / 1 | 0 / 0 / 0 | Low-support S05 candidate；未过 high-confidence gate |
| C2 | 0 / 3 | 0 / 0 | 0 / 0 / 0 | P41 carry-forward candidate；未过 high-confidence gate |
| C3 | 1 / 1 | 0 / 1 | 0 / 0 / 0 | 两个 +doxycycline samples 中均为 low support；未过 gate |
| C4 | 0 / 1 | 0 / 1 | 0 / 0 / 0 | Low-support S05 candidate；未过 gate |

C2 的 S05 raw-support/depth proxy 为 3 / 14.760 = 0.203252；S04 为 0 / 9.627。该比值只作 raw-support/depth proxy，因为 repeat-context evidence 与 standard depth 不能组成经验证的 allele-fraction denominator。S01-S03 的 C2 standard local depth 分别为 21.146x、19.126x、23.101x，均低于预设 60x sensitivity context；其中的 zero raw counts 不能解释为 high-sensitivity exclusion。

### 7.6 客户要求的 genome-wide source-file availability

| Requested data category | Source entries present | 已记录的 index status | Aggregate source size | Availability statement |
| :--- | :---: | :---: | :---: | :--- |
| BAM | 18/18 | 18/18 indexes present | 907,440,957,244 bytes（907.44 GB；845.12 GiB） | Source data 可用；当时不在补充包中复制 |
| Manta filtered SV VCF | 18/18 | 18/18 indexes present | 3,792,358 bytes | Filtered source 可用 |
| Manta candidate SV VCF | 18/18 | 0/18 index entries recorded | 54,016,252 bytes | Candidate-stage source；不改称“unfiltered” |
| Manta candidate small-indel VCF | 18/18 | 18/18 indexes present | 42,142,264 bytes | Candidate small-indel source；不是 project-wide SNV/indel callset |
| DELLY discovery BCF | 18/18 | 0/18 index entries recorded | 30,040,099 bytes | Discovery-stage source 可用 |
| DELLY region-genotyped BCF | 18/18 | 0/18 index entries recorded | 454,914 bytes | Region-genotyped source 可用 |
| Project-wide SNV/indel VCF | 0/1 | 0/1 | 0 bytes | 当前 inventory 中未找到具 provenance 的 source |

## 8. 结论、局限与下一步

### 8.1 Candidate-level conclusion table

| Candidate | 当前状态 | 支持该状态的证据 | 尚未建立的内容 |
| :--- | :---: | :---: | :--- |
| C1 | Low-support S05 candidate | S05 中 1 个 raw 且 1 个 standard-eligible dual-endpoint QNAME | Clone/doxycycline specificity 或 confirmed junction |
| C2 | 优先级最高的 P41 EV0623 carry-forward hypothesis | S05 中 3 个 raw QNAME、5 个 resolved read ends、一致的 2-bp `AT` overlap、in-silico assay design record | Independent breakpoint confirmation、clone specificity 或与 cytogenetic observations 的 causal link |
| C3 | 两个 +doxycycline samples 中的 low-support candidate | 每个 +doxycycline sample 1 个 raw QNAME；S05 有 1 个 standard-eligible QNAME | High-confidence event status 或 sequence-level junction |
| C4 | 伴有 nonconcordant DELLY anchor context 的 low-support S05 candidate | S05 中 1 个 raw 且 1 个 standard-eligible QNAME；重复 DELLY chr11 anchor-only records | chr11-chr3 BND concordance 或 clone/doxycycline specificity |

### 8.2 如何阅读结果：主张与 evidence paths

| Core conclusion | 为什么这样解读 | Evidence location |
| :--- | :---: | :--- |
| C2 获得 base-pair detail，但未获独立 confirmation | 3 个 QNAME 的 5 个 read ends 都为 `AT` 2-bp overlap/0-bp gap；repeat context 与 single WGS data source 构成解释边界 | Section 7.2；`results/c2_read_end_junction_resolution.tsv`；`results/c2_junction_resolution_summary.tsv` |
| 无 exact Manta/DELLY caller concordance | Two-end 50-bp rule 防止将单一 chr11 insertion anchor 混同为 chr11-chr3 BND | Section 7.3；`results/candidate_manta_delly_concordance.tsv` |
| C1-C4 均低于 clone-specific gate | Raw support、standard-eligible support、local depth、sequence-level criteria 分开记录；C2 S05 standard-eligible support 为 0 | Section 5 与 7.5；`results/clone_specific_comparison_with_criteria.tsv` |
| Rescue 对 C1-C4 的增量为测得的 0 | 5 个 recovered unmapped mates 保留在 parent rescue ledger，candidate intersection 明确计为 0 | Section 7.4；`results/read_rescue_candidate_impact.tsv` |
| Transfer planning 可从 manifest 开始，但有一项 gap | Inventory 记录 source path、presence、index、bytes；并未把不存在的 project-wide SNV/indel VCF 说成存在 | Section 7.6；`results/data_delivery_availability_manifest.tsv`；`results/customer_request_coverage.tsv` |

### 8.3 局限

本工作是对同一短读长 WGS alignments 的更深入分析，不是 independent biological replicate 或 assay。Split-read alignment representation 的 2-bp overlap 与 microhomology 一致，但单凭该项不能证明 repeat locus 中唯一的 breakpoint。即便观察到 Manta/DELLY agreement，也只是同一数据的 algorithmic evidence，而非 external confirmation。

当前未提供 independent breakpoint-spanning PCR result、新 sequencing dataset 或 raw fluorescence in situ hybridization（FISH）image/scoring table。现有 local depth 低于 3-fragment/5% criterion 所采用的 60x sensitivity context。因此，parental/control 中的 zero support 不是经验证的 high-sensitivity negative result。当前 inventory 未识别到具 provenance 的 project-wide SNV/indel VCF source。

### 8.4 下一步

对 C2 合理的下一项科学证据是：以合适 control 和实测 analytical sensitivity，在 independent DLD1 genomic DNA 上测试先前选择的 in-silico breakpoint-spanning PCR pair。若另行要求 sensitivity study，应先以预设 60x context 评估实际 local short-read WGS depth，再解释 three-fragment/5% detection target。若要将 genomic candidate 与 cytogenetic subpopulation 联系，必须先提供 raw FISH images 与 per-cell scoring table。

准备本科学报告时没有执行 external data transfer、client contact、commercial quote 或 billing action。客户所述 $840 cap 及任何 transfer method 都须由 Athenomics 另行行政确认，之后才能实施任何对外动作或超出本报告分析范围的工作。

## 9. 补充包文件清单

| 补充包内路径 | 内容 | 用途 |
| :--- | :---: | :--- |
| `P42_Client_Request_Completion_1024.md` | 英文原报告 | 科学审阅入口 |
| `figures/P42_Fig14_Client_Request_Completion.png` 与 `.pdf` | 一张客户问题完成图（两种格式） | 汇总 P41→P42 refinement、rescue impact 与 caller comparison |
| `results/c2_junction_resolution_summary.tsv` | 1 条 C2 summary row | Machine-readable microhomology/gap summary |
| `results/c2_read_end_junction_resolution.tsv` | 5 条 resolved read-end rows | 保留 FLAG、CIGAR、MAPQ、query coordinates 的 alignment-level evidence |
| `results/candidate_manta_delly_concordance.tsv` | 20 条 candidate-by-DLD1-sample comparisons | Exact two-end 与 anchor-only caller classifications |
| `results/clone_specific_comparison_with_criteria.tsv` | 20 条 candidate-by-DLD1-sample rows | Support、depth 与 gate criteria |
| `results/read_rescue_candidate_impact.tsv` | 4 条 C1-C4 rows | Rescue 前/后 accounting |
| `results/customer_request_coverage.tsv` | 5 条 request-to-output rows | Scope coverage 与 interpretation boundaries |
| `results/data_delivery_availability_manifest.tsv` | 109 条 source-inventory rows | 不复制 source data 的 transfer-planning inventory |
| `results/customer_question_summary.tsv` 与 `results/claim_evidence_matrix.tsv` | Compact answer 与 claim-source tables | 快速审阅详细证据 |
| `methods/source_manifest.json` | Source SHA-256 与 inventory summary | Reproducibility 与 source traceability |

---

供科学审阅；未发送给客户。

Zhen Gao, PhD  
Principal Bioinformatics Scientist  
Athenomics
