# P42：开源 short-read WGS SV 工作流比对与独立验证设计

创建日期：2026-10-23  
状态：实施前决策记录；不包含任何 SV 生物学结论。  
适用数据：18 个既有 hg38 BAM；chr11 primary target `114,604,000–114,625,500`，discovery window `114,554,000–114,675,500`。  
关联计划：[P42_research_plan_and_review_1022.md](P42_research_plan_and_review_1022.md)。

## 1. 本次 GitHub 方法调研的结论

调研重点不是寻找一个“替代 Manta/DELLY 的唯一正确软件”，而是检查成熟项目如何分别处理：

1. 全基因组 discovery；
2. 重复区的多重定位、single breakend 和局部组装；
3. 跨样本相同候选的重新分型（regenotyping）；
4. caller 的内部过滤与原始 mapping evidence 之间的边界；
5. 可复核的中间证据和可视化输出。

| 项目 | 实际做法 | 对 P42 的可迁移部分 | 不能照搬的部分 |
| :--- | :---: | :---: | :---: |
| [Manta](https://github.com/Illumina/manta) | 全基因组 breakend graph → 每个候选局部 assembly/评分；discovery scan 明确去掉 secondary、supplementary 与 unmapped records。 | 旧 Manta 结果可作为独立 caller 注释，并可在目标窗口以解除 repeat-edge 上限的配置重跑。 | 被 Manta 跳过的 record 不是“没有 SV 价值”；不能把它的内部输入过滤变成 P42 原始证据的过滤规则。 |
| [GRIDSS2](https://github.com/PapenfussLab/gridss) | alignment-guided positional de Bruijn graph assembly，支持 single breakend；靶向模式先提取 ROI overlap fragment 的 mate 与 supplementary records，再组装。输出保留 multi-mapping unique allocation 前的 raw breakpoint VCF。 | 最贴合重复区：并列保存 raw breakpoint、unique-allocated support、single-breakend 和 assembly sequence，并进行 RepeatMasker 注释。 | GRIDSS 要求每个 read 只有一条主要 alignment（supplementary 可以，secondary 不可以）；因此只能使用单独生成的 caller-compatible 输入，绝不能取代原始全 FLAG 审计轨。 |
| [SvABA](https://github.com/walaj/svaba) | 可在 BED/区域内做局部 assembly，BWA-MEM 重比对 contig，输出 VCF 及含每样本完整证据的 `bps.txt.gz`。 | 是对本项目 MEGAHIT + all-placement contig alignment 的独立算法复核；适合对同一目标区逐样本产生可比断点评分。 | 其 caller score 不可等同于真实位置；重复区 alternate placement 仍须保留、报告并判为歧义。 |
| [UCLA call-gSV](https://github.com/uclahs-cds/pipeline-call-gSV) | DELLY/Manta + CNV discovery，随后对已发现 SV/CNV 进行 regenotyping。 | clone-specific 判定必须在固定候选坐标上逐样本重测，不能以“某样本 VCF 没有一行 call”代替阴性证据。 | 原项目主要是 germline/CNV 管线，不能直接给 DLD1 clone 特异性设定阈值。 |
| [nf-core/sarek](https://github.com/nf-core/sarek) | 容器化、输入 QC、多个 caller、VCF 标准化/注释和 QC 汇总；可做共识。 | 可复现的 manifest、版本、调用参数、per-sample QC 与标准化 caller 字段应成为交付的一部分。 | 多 caller 交集不能作为重复区的硬性淘汰条件，否则会把真但低 mappability 的 junction 误删。 |

上表的 Manta 过滤行为来自其[方法文档](https://github.com/Illumina/manta/blob/master/docs/methods/primary/methods.tex)；GRIDSS 的 targeted extraction、one-alignment 输入要求、raw-before-allocation VCF、single-breakend 定义和 joint calling 建议来自其[官方 README](https://github.com/PapenfussLab/gridss)及[QuickStart](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md)；SvABA 的 targeted run、contig realignment 和 verbose evidence 来自其[官方 README](https://github.com/walaj/svaba)。

## 2. 与当前 P42 设计直接相关的关键区别

### 2.1 “caller 输入”与“原始证据 universe”不能混为一谈

Manta 的发现阶段排除 secondary/supplementary/unmapped alignment 是其计算模型的一部分，不是关于这些 records 生物学价值的判断。在本项目的重复区问题中：

- **secondary** 是同一 read 的 alternative full alignment；它不能增加独立分子数，却是候选真实位置不唯一的直接证据；
- **supplementary** 是同一 read 的 split/chimeric alignment；它不能增加独立分子数，却可能直接跨越 junction；
- **unmapped mate** 可能保存一端无法直接定位、但另一端位于目标区的实际序列；
- MAPQ0/multi-mapping 也不能自动视为无关，因为研究问题正是“在相似位置之间哪一个是真实来源”。

因此，P42 的 full-evidence track 继续保留并审计所有这些类别。GRIDSS/Manta/SvABA 的 caller-compatible track 是**第二条平行证据轨**，其“每 read 一条代表 alignment”规则只用于让该 caller 合法运行与得到可解释的 caller score。

### 2.2 “unique allocation”是计数规则，不是删除原始位置

GRIDSS 明确把 raw breakpoint call 与 multi-mapping read 的 unique allocation 分开。P42 采用同一思想：

- 原始 ledger 永远保存 `QNAME`、read-end、FLAG、reference/start/end、CIGAR、MAPQ、SA/其他定位及进入/未进入某一 caller-compatible 输入的原因；
- 主支持数为独立 `QNAME/template` 数；
- 并列报告 `QNAME + read-end` 的 read-level 数；
- 如果一个 template 有多处 mapping，它只在 template count 中计一次，但其所有 mapping locations 仍保留作定位歧义审计；
- caller 的 allocated read support 只能解释为“该算法在候选之间的归属”，不能解释为“其他 mapping 已被证伪”。

### 2.3 caller 缺失不构成阴性

已有 Manta/DELLY VCF 的 “not called” 只能描述为 *not called under that caller/version/configuration/filtering rule*。特别是 Manta 的默认 `graphNodeMaxEdgeCount=10` 会在两端都高连通的 repeat graph edge 上静默跳过候选；P41 的历史 Manta 配置将单独审计。若未确认使用提高上限的配置，则以新目录运行靶向、repeat-unmasked Manta 作为补充，而不覆盖既有全基因组 VCF。

## 3. P42 新增并行验证轨的 decision table

| 研究问题 | data universe / 观察单位 | 保留与去重规则 | 判定指标 | 允许的结论 |
| :--- | :---: | :---: | :---: | :---: |
| 在目标窗口发现哪些异常 mapping/junction 线索？ | 18 样本目标/发现窗口内的全部 records，以及 target-overlap fragment 的 mate；record、read-end、template 三层均保留。 | 不用旧 caller/read-name list 限定 discovery；所有 FLAG 类别留在 ledger。template 主计数按 QNAME，read-level 按 QNAME+end。 | FLAG、MAPQ、CIGAR/soft clip、SA、mate 状态、pair geometry、定位集合、原始与过滤后数量。 | “观察到/未观察到某类证据”；不能由此单独宣布唯一断点或 clone 特异性。 |
| 某一候选是否能被短读局部 assembly 支持？ | 由 full-evidence track 产生的 59 个 locus；MEGAHIT FASTQ 仅为组装输入，完整 BAM ledger 另存。 | FASTQ 每个 QNAME/read-end 选择最长的可用 stored SEQ；不删除其余 alignment ledger。5 个索引无法直取的 unmapped mate 已全 BAM 顺序扫描回收。 | contig 是否含双侧 sequence；all-placement mapping、primary/secondary/supplementary placement、microhomology/insert sequence 和 source-template 支持。 | “assembly-supported / unresolved / repeat-ambiguous”；primary contig placement 不等于真实 locus。 |
| 同一候选在独立 caller 中是否重现？ | 目标窗口的 caller-compatible 衍生输入；DLD1 S01–S05 是 clone-specific 主比较，18 样本均保留为背景/可见性上下文。 | 为满足 caller 要求制作新输入；secondary 不喂给 GRIDSS，但在原始 ledger 保留并输出排除原因。supplementary/mate 依工具规范处理。 | GRIDSS raw vs allocated call、single breakend、assembly score、SvABA VCF/bps evidence、unmasked-Manta call；每项标记运行状态。 | “algorithmically concordant / discordant / not evaluable under a caller”；不是独立生物学验证。 |
| 某事件是否 Dox clone-specific？ | 固定 breakpoint/candidate coordinate 在 S01–S05 的逐样本 re-evaluation；其余 13 样本仅作跨细胞系背景。 | 候选先由全证据 + assembly/caller union 生成，再在所有 5 个样本按相同窗口与阈值重新计数；不以 VCF 缺失代替 0。 | 每样本 template/read-end support、reference support、coverage、MAPQ/placement ambiguity、caller/assembly status、S04/S05 与 S01/S02/S03 的差异。 | “Dox-enriched/clone-specific candidate”“shared”“insufficient evidence”“ambiguous”；不把 read rescue/re-alignment称为独立实验验证。 |

## 4. 已完成的输入与组装完整性检查

截至本文件创建：

- 全队列 raw evidence、caller audit、coverage 与 corrected QNAME/template clustering 已完成；`51,364` 是 evidence bin 数，不是 51,364 个 SV；
- 59 个局部 assembly locus 共纳入 2,004 个 priority QNAME；
- 5 个 `bam.mate` 显示 “mate is unmapped” 的 cases 已扫描相关完整 coordinate-sorted BAM 至 EOF，均找回 1 个 R2、各长 150 nt 的 unmapped sequence；
- 59/59 locus MEGAHIT 完成，0 failed、0 insufficient-sequence；组装生成 90 条 contig、合计 121,426 bp；
- 所有 contig 已针对 hg38 使用 all-placement alignment；该结果必须先做 multi-placement/ambiguity 分类，尚未据此输出任何 SV 结论。

相关 machine-readable 审计文件：

- `analysis/20261023_fullcohort_repeat_resolve/05_local_assembly/inputs_v1/unmapped_mate_rescue/unmapped_mate_recovery_status.tsv`
- `analysis/20261023_fullcohort_repeat_resolve/05_local_assembly/assemblies_v1/audit.json`
- `analysis/20261023_fullcohort_repeat_resolve/05_local_assembly/assemblies_v1/assembly_run_status.tsv`
- `analysis/20261023_fullcohort_repeat_resolve/05_local_assembly/contig_placements_v1/contig_placement_summary.tsv`

## 5. 将实施的新增步骤（均为新目录，不覆盖历史产物）

1. 建立独立环境 `p42_sv_validation_20261023`，安装 GRIDSS2 与 SvABA；记录精确版本和二进制自检。
2. 审计 P41 的 Manta 配置，明确历史全基因组 caller 是否采用默认 repeat-edge 限制；原文件只读保留。
3. 写出可复现的 numbered scripts，生成 caller-compatible targeted inputs，并为每个未纳入该输入的 record 输出 filter-impact audit。完整 full-evidence BAM ledger 不变。
4. 对 fixed target / discovery regions 运行：repeat-unmasked targeted Manta（仅在历史配置未证实已解除限制时）、GRIDSS targeted joint/per-sample validation、SvABA targeted local assembly。每个工具独立目录、独立日志、独立 `.complete` 标记。
5. 标准化候选及两个断点端，建立 raw-support、unique-allocated support、contig support、caller status 和每样本 regenotyping 矩阵。
6. 对每个高优先级事件制作多证据图：coverage、template-level read geometry、FLAG/mapping class、contig all-placement、caller concordance、DLD1 five-sample comparison。图中必须显式标记 mapping ambiguity，不把同一 template 的多条 alignment 画成独立分子。

## 6. 预先规定的限制表述

1. 反复对同一 short-read BAM 重提取、局部组装和重调用，只能检验 mapping/algorithmic reproducibility；不是新的独立实验或独立 DNA 分子。
2. 在强重复区，能够给出“多处候选位置集合”和“某位置更受支持”的证据等级，但不承诺短读一定可唯一定位每一端。
3. caller concordance 提升可信度但不是必要条件；caller discordance 也不是自动否定。二者都必须回到 template-level 原始证据、contig placement 和 control comparison 解释。
4. 没有可用 sequence、没有足够独立 template、或 contig 只能多处等价比对时，结论必须为“unresolved/ambiguous”，而不是“absence of SV”。
