# P42 当前 Analysis 1：已发现问题与重跑前审计

- **文档性质**：内部技术审计与重跑设计输入。记录当前 Analysis 1 已发现的方法学不足，不替代任何结果表，不删除或改写已生成的历史产物。
- **创建日期**：2026-10-22
- **当前状态**：用户已完成审阅并授权按严格新设计分析全 18 个样本。历史 candidate-gated 结果已归档保留；新版 full-window inventory、caller audit、coverage 与 18 样本 raw junction-evidence extraction 已完成。本文件继续记录新旧主线切换中发现的任何 evidence-model 问题；在本文件标为“待确认”的实质口径解决前，不将新版中间产物表述为完整答案。
- **适用范围**：全 18 个样本、GRCh38 `chr11:114,604,000–114,625,500` target window、`chr11:114,554,000–114,675,500` discovery window 及其 read-level evidence。DLD1 主比较为 S04/S05 对 S01/S02/S03。

## 1. 审计结论

已完成的 P42 Analysis 1 是一个**限定口径的 candidate-supported-read rescue**：它从 P41 已保留的 179 个候选出发，重新提取这些候选的 source sequence，并以 `bwa-mem2 -a` 重比对后测试固定的两个 endpoint 是否再次出现。它可以回答“这 179 个预选候选中的原 sequence 能否在给定阈值下再次形成记录的 endpoint pattern”。

它**不能充分回答**本项目更广的研究目的：在 repeat-affected / ambiguous chr11 区域，系统审计窗口内所有 read 的竞争定位，并据此判断哪些位置具有更强的真实定位证据。原因不是 `bwa-mem2` 或 BAM 数据本身失效，而是旧流程在 candidate universe、输入 record 保留和 alternative-mapping 解释三个层面作了过窄的选择。

因此，现有结果只能作为历史性的、条件化的 read-rescue 记录保留，不能被作为重复区真实位置判定或完整 junction discovery 的最终答案。该结论不删除或否定表中已发生的重比对事实；它明确这些事实所能支持的推断边界。

## 2. 已执行流程实际回答的问题

| 分析环节 | 已执行的规则 | 实际得到的证据 | 不能据此推出的结论 |
| :--- | :---: | :---: | :---: |
| Candidate universe | 仅保留 P41 `in_target == True` 且 `n_split_nondup_mapq20 >=1` 的 junction | 179 个已知候选的列表 | 窗口内不存在其他 junction、或 P41 未保留 read 没有信息 |
| Read retrieval | 仅查询 P41 `source_reads` 中列出的 QNAME | 这些 QNAME 的可用 `SEQ` 被写入 FASTQ | 窗口内其他 read 与 junction 无关 |
| BAM record filtering | 排除 secondary (`0x100`)；保留常规 supplementary (`0x800`)；同一 `(QNAME, mate)` 仅选一个 sequence record | 每个候选取得一个可重比对的 source-sequence representation | 原始 BAM 的竞争 placement landscape 已完整保存或被系统比较 |
| Independent remapping | `bwa-mem2 mem -a -Y` 保留新的 alternative alignments | 被选取 sequence 的可重比对性及固定 endpoint 是否匹配 | 每个 alternative placement 的相对可信度、真实 locus 或完整 SV architecture 已确定 |
| Candidate classification | 两端 endpoint 在 ±20 bp 内匹配，且两段 `MAPQ >=20` | 4 multi-fragment、173 single-fragment、2 not reconstructed | 单 fragment 即为 confirmed SV，或 primary placement 就是真实位置 |

## 3. 已发现的问题

### 3.1 问题 A：candidate universe 先被 P41 结果封闭

旧流程的第一步并不是对窗口内的原始 read 做 discovery，而是从 P41 target-evidence table 中挑出 179 个已有 junction。凡是不满足 P41 `in_target` 或既定 split-read/MAPQ 条件的 read，都不会进入候选表，更不会进入后续 sequence extraction 与 remapping。

在普通高唯一性区域，这种“从现有候选做 validation”的设计可以合理地节约计算并回答 narrow validation question；但本项目目标区已知为 repeat-affected / ambiguous region，研究问题恰恰包括哪些 mapping placement 才更可信。因此，P41 没有保留的 window read 本身就是必须审计的潜在证据，不能预先视为无关。

**影响**：现有 179 的分母不是“窗口内所有证据”的分母，而是“P41 已保留且满足旧 eligibility rule 的候选”的分母。它不提供 discovery sensitivity，也不提供未入选 read 的审计记录。

### 3.2 问题 B：`read.query_name not in wanted` 直接排除了窗口内绝大多数 read

旧提取脚本的条件为：

```python
if read.query_name not in wanted or read.is_secondary or not read.query_sequence:
    continue
```

其中 `wanted` 只含 P41 已保留 candidate 的 `source_reads` QNAME。一个 record 即使实际落在查询窗口内，只要其 QNAME 不在 `wanted`，就不会被记录、计数、导出 FASTQ 或检查其 partner placement、CIGAR、clip、MAPQ 与 junction 关系。

这与“窗口内全部 read 值得重点研究”的研究原则不相符。窗口坐标是本轮空间筛选条件；是否已被旧 candidate table 引用不能成为丢弃原始 evidence 的理由。

**影响**：旧流程没有留下“窗口内总 record 数、总 unique QNAME 数、各 FLAG 类别数、未进入 P41 candidate 的 read 数及其定位特征”的可审计底表。因而不能量化旧候选筛选遗漏了什么。

### 3.3 问题 C：secondary alignment 被当作可直接丢弃的 record，而非需要保留的定位证据

`read.is_secondary` 是 FLAG `0x100`。它表示同一 read 存在竞争性 alternative mapping；在重复序列区域，这种竞争位置正是判断“该 read 是否能被可靠定位到当前 locus”的关键信息。它不是新的独立 DNA fragment，不能重复计入 molecule support；但“不能重复计数”与“可以不保存、不审计”是两件不同的事。

旧流程的合理局部动机是：同一个 QNAME 的 bases 只需要导出一次用于 `bwa-mem2` 重比对，并且后续 `-a` 会输出新的 alternative alignments。这个动机只足以避免重复导出同一 sequence，**不足以**支撑在原 BAM evidence 层面删除 secondary placement。若某个 read 在窗口内仅以 secondary record 形式出现，旧代码会直接跳过它；若 secondary 与 primary/supplementary 同时存在，旧结果也不会保存原 BAM 中它们的坐标、MAPQ、CIGAR 及彼此竞争关系。

**影响**：现有分析不能系统比较“primary location 与 secondary location 哪个更支持目标 junction”，也不能识别“read 在多个相似位置间本质上不可判定”的情形。对真实 locus 的推断需要保留全部 alignment evidence，再在 QNAME/mate 层面去重，而不是提前删除 FLAG 类别。

### 3.4 问题 D：supplementary 与 secondary 的证据角色没有在输出中被分开审计

supplementary alignment（FLAG `0x800`）通常表示同一 read 的另一非重叠 segment，是 split-read / chimeric structure 的组成部分；secondary alignment（FLAG `0x100`）通常表示同一 read 的竞争 placement。二者都不应当被当作额外独立分子，但它们回答不同的科学问题：supplementary 可描述 junction segment；secondary 可描述 repeat-induced localization ambiguity。

旧脚本保留常规 supplementary record 的 sequence，但没有生成一个 raw alignment audit table，逐 read 显示 primary / supplementary / secondary 分类、FLAG、RNAME、POS、MAPQ、CIGAR、SA tag、clip 与 candidate endpoint 的关系。因此，即使部分必要 record 曾被读到，也无法从现有交付物复核“为什么选择该位置为更可信定位”。

**影响**：输出将“可再次匹配固定 endpoint”压缩为 candidate status，丢失了支持该 status 的完整 competing-alignment context。

### 3.5 问题 E：`-a` 的使用是必要但不充分的

`bwa-mem2 mem -a -Y` 的 `-a` 确实会输出重新比对时找到的 alternative alignments；这是旧流程中正确的防止“只显示一个 best hit”措施。但随后分类脚本只寻找能匹配预设两端 endpoint 的 alignment pair，并在多个 pair 中选择 alignment score (`AS`) 之和最高的一对，再将结果压缩为 unique / ambiguous / not reconstructed。

这种设计没有输出或比较每个 competing locus 的完整集合、score differential、MAPQ distribution、原 BAM 与再比对结果的一致性、或多个候选 endpoint 之间的竞争关系。`MAPQ >=20` 是一个 operational threshold，不是重复区域真实位置的证明。

**影响**：保留 `-a` 不自动等于完成 ambiguity adjudication。新的流程必须将所有 relevant alignment 显式制表和解释，而非仅把最终符合预设 endpoint 的 pair 留在分类结果中。

### 3.6 问题 F：`SEQ` 非空不等于完整原始 read 已恢复

旧脚本只要求 `read.query_sequence` 非空。该条件正确地区分了“有可重比对输入”与“完全没有 stored bases”，但不会判断该 record 是否含完整原始 read sequence。带 hard clipping（CIGAR `H`）的 record 仍可能有非空 `SEQ`，但被 hard-clipped 部分 bases 不在该 record 中。脚本优先选择窗口内 primary record；若 primary 不在窗口，则可选用 supplementary record。这对有限 read-rescue 可用，但不足以保证每一个输入都适合 base-pair-level breakpoint interpretation。

**影响**：下一版必须逐 record 审计 CIGAR、hard/soft clipping、`SEQ` length 与 read-end identity；若需要完整 read，应回溯原始 FASTQ 或具有完整 sequence 的 record，而不能仅以 `SEQ` 非空作为充分条件。

### 3.7 问题 G：过滤影响缺少逐类量化审计

旧流程的最终 QC 证明 179 个预选 candidate 全部有至少一个可评估 source fragment，但没有逐类输出：窗口内总 record 数、非目标 QNAME 被排除数、secondary 被排除数、缺 `SEQ` 被排除数、primary/supplementary 选择/替换数，以及这些集合在 sample 和 read-name 层面的重叠关系。

**影响**：没有这份 audit，无法区分“某类 evidence 根本不存在”与“它存在但被代码过滤”。严谨的过滤不仅要写条件，还要报告每个条件对 data universe 的实际影响。

### 3.8 问题 H：read-end 去重不能替代 fragment/template 去重

新版预注册表最初把 `QNAME+read-end` 称为 molecule。这一规则对 **alignment context 去重**是必要的：同一 read-end 的 primary、secondary、supplementary records 不得被重复计成支持。然而 paired-end WGS 中，同一 QNAME 的 R1、R2 通常来自同一个 library template/fragment；两端都落入同一 cluster 时，若把它们直接加为两个“independent molecules”，会在 fragment-level support 和 AF 阈值上放大证据。

这不意味着必须丢弃任一 mate：两端可提供互补的 orientation、insert size 和跨 breakpoint 信息，仍须完整保留。正确做法是并列两个不可混用的量：

- **independent read-end count**：`(sample, QNAME, read-end)` 去重，回答有多少不同 sequenced read ends 给出该证据；
- **independent fragment/template count**：`(sample, QNAME)` 去重，回答有多少独立 library fragments 支持该 cluster；高置信 SV 的“≥3 independent molecules”应以此口径评估。

**处理与当前状态**：cluster 输出同时写两列、逐 row 保留 QNAME/read-end/FLAG。用户于 2026-10-23 确认：`(sample,QNAME)` fragment/template count 是 `≥3 independent molecules` 的主阈值，`(sample,QNAME,read-end)` 作为 read-level 支持数并列报告；据此归档未运行的旧脚本、修订脚本并恢复 cluster 运行。此前暂停是记录明确的科学口径阻塞，不是无人监控或无理由空转。

## 4. 对现有结果应保持的正确表述

下列已运行的事实仍然可以保留，但必须带上完整限定语：

| 旧结果表述 | 可保留的严格表述 | 不可扩展为 |
| :--- | :---: | :---: |
| 177/179 candidates reconstructed | 在 P41 预选的 179 个候选、被旧脚本成功选取并重新比对的 source sequence、既定 reference 和 ±20 bp / `MAPQ >=20` 判据下，177 个候选有至少一个 endpoint-pair reconstruction | 窗口内 177 个真实 SV、或该区域只有这 179 个候选 |
| 4 multi-fragment supported | 在旧 candidate universe 内，有 4 个候选被至少 2 个不同 QNAME 的 selected source sequence 重建 | 4 个已确认的重排或完整断点 |
| 173 single-fragment supported | 在旧 candidate universe 内，有 173 个候选仅有 1 个不同 QNAME 的 selected source sequence 重建 | 173 个可直接用于克隆、处理效应或生物学差异推断的事件 |
| 2 not reconstructed | 在旧固定 endpoint 和 selected source-sequence 输入下，这 2 个候选未形成合格 endpoint pair | 两个事件在生物学上不存在、或 P41 记录必然错误 |

## 5. 重跑前必须完成的新设计

用户尚未授权重跑。待全文审阅完成并明确启动后，新的 workflow 在写任何正式脚本和执行前必须具有以下可审查设计：

1. **Window-wide record inventory**：对每个 S01–S18 BAM 读取目标窗口及经定义的必要 flank 内的全部 primary、secondary、supplementary record，而不是从 P41 candidate/read-name 列表起步。输出 raw audit table，至少包含 sample、QNAME、mate、FLAG / primary-secondary-supplementary class、RNAME、POS、MAPQ、CIGAR、`SEQ` length、hard/soft clip、SA tag 和 duplicate/QC state。
2. **Evidence、read-end 与 fragment/template 分层**：保留每一条 alignment record 作为定位 evidence；以 `(sample, QNAME, read-end)` 合并同一 read-end 的 alternative records，以 `(sample, QNAME)` 统计独立 fragment/template support。绝不以“避免重复计数”为由删除 alignment context。
3. **Candidate discovery 不依赖旧候选 gate**：从 window-wide evidence 建立新的 junction / competing-locus candidates；P41 179 个候选只能作为 comparison annotation，不可作为唯一输入。
4. **Repeat-aware localization**：对同一 read 的 primary、secondary、supplementary placement 显式比较，而不是将 primary label或单一 `MAPQ` 当作真实 locus。报告所有竞争位置、相对 alignment score、MAPQ、CIGAR / clipping、endpoint overlap 及其对定位结论的影响。
5. **Sequence completeness audit**：识别 `SEQ='*'`、hard clipping、缺 quality 和 read-end 缺失；对于需要完整 sequence 的 read，定义回溯原始 FASTQ 或其他完整 record 的规则，并报告未能恢复的数量。
6. **Pre-registered decision table**：在运行前写清每类 evidence 的 input definition、过滤条件、保留原因、计数单位、输出字段、可以支持的结论及不能支持的结论；用户可在正式运行前审阅该表。
7. **Filter-impact QC**：每一个过滤条件必须输出 before/after 数目、被排除 read/QNAME 的可追踪清单或可复算汇总，不能只输出通过后的最终集合。
8. **独立验证边界**：明确区分“同一 source sequence 的重比对复现”“多分子支持”“原 BAM mapping evidence”“可判定真实 locus”与“需要 orthogonal confirmation”的结论层级。

## 6. 历史产物保护与主线切换原则

当前结果、脚本、报告和运行日志均保留，不删除、不原地覆盖。待新设计获审阅并完成重跑后：

- 旧 candidate-supported-read rescue 的脚本、analysis outputs 和内部说明移入项目内 `superseded_for_comparison/` 的独立版本目录，并附 README 说明旧问题、运行口径和与新版的差异；
- 已生成的 supplementary delivery 作为历史记录保留，不能被新结果覆盖；
- 新 workflow 使用新的带日期目录、新的 numbered scripts、独立的 machine-readable audit tables 和新的 report source；
- 新旧结果必须并存，报告不得将旧结果写成被静默删除或不可追溯的“错误版本”。

在用户明确“全文审阅完成，可以重跑”之前，上述内容只是设计约束，不执行归档移动或新分析。

## 7. 关联文件

| 文件（相对于项目根目录） | 本审计中的角色 |
| :--- | :---: |
| `docs/P42/P42_chr11_junction_read_rescue_method_and_conclusions_1022.md` | 旧流程的详细方法与条件化结论说明；已补充 FLAG / `SEQ` 边界 |
| `docs/P42/P42_research_plan_and_review_1022.md` | 当前分析指南、计划与审阅真源；其中保留旧流程的问题审计与重跑判据 |
| `scripts/20261022/1_prepare_candidates.py` | 旧 candidate universe 的构建逻辑 |
| `scripts/20261022/2_extract_support_reads.py` | 旧 read retrieval、secondary exclusion 和 sequence selection 逻辑 |
| `scripts/20261022/3_classify_rescue.py` | 旧 endpoint-pair classification 逻辑 |
| `analysis/20261022/results/rescue_candidate_results.tsv` | 旧 179-candidate result；仅能按第 4 节限定语解读 |
