# P42 BAM alignment 与 junction evidence 术语说明

**用途：** 本文定义 P42 在 chr11 repeat-aware SV/CNV 分析中使用的 BAM/SAM alignment 与 junction-evidence 术语，帮助阅读 alignment ledger、候选表、IGV 截图和后续图件。

**适用范围：** 这是内部判读辅助材料，不是新的分析规则或生物学结论。P42 的证据 universe、计数口径、过滤审计及可允许结论以 [`P42_research_plan_and_review_1022.md`](P42_research_plan_and_review_1022.md) 为唯一真源。

## 1. 一条 BAM alignment record 包含什么

一条 BAM record 是一条测序 read 的**一次比对记录**，而不是一个独立分子，也不自动等同于一个真实基因组位置。P42 重点读取以下字段：

- **QNAME**：read 名称。同一 paired-end library fragment 的 R1 和 R2 通常共享 QNAME。
- **FLAG**：位标志；例如 `0x100` 是 secondary、`0x800` 是 supplementary。一个 FLAG 可同时含多项状态。
- **RNAME / POS / strand**：此记录所报的参考 contig、位置和方向。
- **CIGAR**：read 序列如何与该位置对齐，例如 `70S80M`。
- **MAPQ**：该 alignment 的 mapping quality；低值或 0 常提示重复区/多重定位歧义，不能自动视为无价值或错误。
- **SEQ**：BAM 中保存的 read 序列。`SEQ` 非空仅表示可用于重比对或局部组装；不能自动证明它是完整原始 read。
- **SA:Z tag**：同一 read 的其他 supplementary alignment 摘要；见第 3 节。

## 2. CIGAR、soft clip 与 hard clip

常见 CIGAR 操作含义如下：

- **`M`**：read 与参考的对齐碱基（可包含 match 或 mismatch）。
- **`I` / `D`**：相对参考的 insertion / deletion。
- **`S`（soft clip）**：该段 read 序列未纳入当前 alignment，但仍保存在 `SEQ` 中。例如 `70S80M` 表示前 70 nt soft-clipped、后 80 nt 在当前位置比对。
- **`H`（hard clip）**：该段序列不在此 BAM record 的 `SEQ` 中。因此不能由这条 record 单独恢复其碱基序列。

soft clip 是“当前位置没有解释 read 的这一段序列”的现象，不是 SV 的同义词。它可由真正 junction、低复杂度、重复序列、测序末端错误或局部比对失败造成。P42 中，soft-clipped sequence 只有在其替代定位、方向、独立 fragment 支持和（如有）assembly 结果一致时，才可能提升为 junction evidence。

## 3. SA、split read 与 supplementary alignment

**split read** 指同一条 read 的不同连续片段分别比对到不同位置或不同方向；这常由一个结构变异断点、融合或嵌合结构引起，也可由重复区歧义引起。

**supplementary alignment** 是该 split read 的补充比对记录，FLAG 含 `0x800`。通常同一 read 有一条 primary record，其余互补片段为 supplementary records。

**`SA:Z`** 是 SAM/BAM optional tag，列出该 read 的其他 supplementary alignment 摘要。其单个条目形如：

```text
SA:Z:chr11,114600000,+,70S80M,60,2;
```

依次表示：`contig, position, strand, CIGAR, MAPQ, mismatch count`。这条信息说“同一 read 还存在一个 alignment 片段”，并不单独证明该两个位置之间存在真实 DNA junction；应回查所有同 QNAME/read-end 的 records、完整 CIGAR/SEQ、替代 placement 与独立 fragment 支持。

## 4. primary、secondary 与 supplementary 的区别

- **primary alignment**：aligner 为该 read 标为主的记录。它是默认代表记录，**不是**真实 locus 的充分证明。
- **secondary alignment（FLAG `0x100`）**：同一整条 read 的替代比对位置，多见于重复区或同源序列的 multi-mapping。它不代表额外独立分子，但可能正是定位歧义的直接证据。
- **supplementary alignment（FLAG `0x800`）**：同一 read 的互补片段的 split/chimeric alignment。它不代表额外独立分子，但可能保留跨 junction 的结构信息。

P42 的 full-evidence track 保留这三类 record 并分层报告。只有为满足 GRIDSS 输入要求而新建的 caller-compatible 派生 BAM 才排除 secondary；原始 alignment context 不删除，并有 filter-impact audit。

## 5. paired-end、discordant mate 与 unmapped mate

- **paired-end fragment/template**：一对 R1/R2 通常来自同一个文库分子。P42 的独立 molecule 主计数单位是 `(sample, QNAME)`，不能将 R1、R2 或同一 read 的多个 records 加总为多个分子。
- **read-end**：R1 或 R2 的单端单位。`(sample, QNAME, read-end)` 用来合并同一 end 的 primary/secondary/supplementary records，作为辅助 read-level 指标。
- **discordant mate/pair**：R1、R2 的距离、方向、contig 组合或 insert size 与文库正常分布不符。它可提示 SV，但也可能来自重复区错配或建库伪影。
- **unmapped mate**：一端有可靠 anchor，另一端未比对。未比对的一端可携带 junction 或外源序列；它不是“没有证据”。P42 会在可用 BAM 记录中进行 mate rescue，并对无法恢复的情况单独审计。

## 6. MAPQ、multi-mapping 与重复区解释

**MAPQ** 是 aligner 对“此 alignment 相对于竞争位置”的质量估计，不是测序碱基质量，也不是 SV 置信度。`MAPQ=0` 常见于多个同等或近似同等的 placement。

在 P42 的重复区问题中：

- MAPQ0、secondary 或多重 placement 不可单独用于声明唯一断点或 copy number；
- 但它们不能因为“不唯一”而从原始 evidence universe 消失；
- 需要同时展示替代位置、方向、CIGAR、SA/split context、独立 template 支持和局部 assembly 的 all-placement 审计；
- coverage 亦需分 standard 与 inclusive track，避免将 mappability 偏差误写成真实 CN 改变。

## 7. P42 的最低判读边界

- 单端 soft clip、单条 SA、caller 的单一记录或 primary placement，只能构成 candidate/ambiguous evidence。
- 同一 QNAME 的 primary、secondary、supplementary records 是同一测序分子的不同 mapping context，计数时必须归并而不是相加。
- 重比对、局部 assembly 或不同 SV caller 的一致，支持算法/定位复现；它们不是独立生物学验证。
- “未观察到”须连同可评估 coverage、有效分子数与检测下限报告；不能直接写成事件不存在。

## 8. 阅读 P42 输出时的最短路径

1. 先看 full-evidence ledger：同一 QNAME/read-end 是否存在 primary、secondary、supplementary、SA 或 mate context。
2. 再看 fragment/template 计数：支持是否来自多个独立 QNAME，而不是同一 fragment 的重复 records。
3. 检查替代 placement、MAPQ 与 CIGAR：是否存在重复区等价解释。
4. 最后整合 caller-compatible 结果、contig all-placement 与 DLD1 对照比较，按当前研究计划的阈值与结论等级报告。

## 9. 相关 P42 产物

- `analysis/20261023_fullcohort_repeat_resolve/` 下的 full-window inventory、junction-evidence cluster、local-assembly ledger 和 contig all-placement 记录保存原始或可追溯的 evidence context。
- `analysis/20261023_fullcohort_repeat_resolve/06_external_sv_validation_v3/` 是正在生成的外部 caller-compatible 轨；其 raw targeted BAM、secondary-exclusion audit 与后续 GRIDSS/SvABA 输出须与 full-evidence track 并列解释。

## 10. 修订记录

- 2026-10-23 — 首次建立；响应对 BAM soft clip、SA/split、supplementary 与 secondary alignment 的术语解释需求。本文不改变既有分析脚本、候选定义或结论阈值。
