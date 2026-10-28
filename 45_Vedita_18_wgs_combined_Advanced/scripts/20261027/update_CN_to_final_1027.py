#!/usr/bin/env python3
# Bring docs/P45_..._1027_CN.md (translated from English c32b1d8) in line with the final 1027 English report
# (user: 审阅完毕, 2026-10-27). Mirrors the 7 English hunks since c32b1d8: header/signature hard breaks, Section 2 /
# Section 8.2 swap with the strengthened list, S05 split into three rows, Methods table names, Section 7.4 three readings.
import re, sys

CN = sys.argv[1]
t = open(CN).read()

# 0. preamble: record which English version this mirrors
t = t.replace("（2026-10-25 起；2026-10-27 同步至 1027 版）",
              "（2026-10-25 起；2026-10-27 同步至 1027 最终版，即用户“审阅完毕”的版本）", 1)

# 1. hard line breaks: header field lines and the signature name line
L = t.split("\n")
for i in range(len(L) - 1):
    if (i < 20 and L[i].startswith("**") and L[i + 1].startswith("**") and ":**" in L[i + 1]) or L[i] == "**Zhen Gao, PhD**":
        assert not L[i].endswith("\\")
        L[i] += "\\"
t = "\n".join(L)

# 2. Section 2 / Section 8.2 swap
s2 = t.index("## 2. 核心结论（Core Conclusions）\n"); s3 = t.index("## 3. 关键发现")
p2 = t.index("### 8.2 实用解读（Practical reading）\n"); p3 = t.index("### 8.3 局限")
old_s2 = t[s2:s3]
m = re.search(r"\n\n(\| 样本组 \| 样本 \|.*?)\n\n---\n", old_s2, re.S)
table_notes = m.group(1)

old_s05 = next(l for l in table_notes.split("\n") if l.startswith("| **DLD1 +Dox 克隆 4** | S05 |"))
new_s05 = "\n".join([
    "| **DLD1 +Dox 克隆 4**：总体 | S05 | 断崖以远丢失的是与 S04 相同的那条 homolog（第 7.7 节）；近端呈渐变的等位基因失衡。断崖附近有两个 junction，S01–S03 中均不存在（见下两行）。 | 断端比 S04 更复杂；这些特征的组合符合在该位点发生反复断裂与重连。 |",
    "| **DLD1 +Dox 克隆 4**：fold-back inversion | S05 | 位于深度断崖内的 376 bp fold-back inversion：chr11:114,605,973 与 114,606,349 头对头相连。 | 这 376 bp 不是插入的外来序列：chr11 读到 114,605,973 后，从 114,606,349 开始反向继续，即断端折回到自身。因此该 homolog 在此点以远没有 chr11 序列，与远端丢失吻合。 |",
    "| **DLD1 +Dox 克隆 4**：chr21 junction | S05 | chr11:114,607,192 到 chr21:20,648,241 的 junction，位于断崖远端 1.2 kb；是四个 caller 中有三个调出的那个事件。 | 连接到 chr21 的只有截至 114,607,192 的 chr11 序列；其余约 20 Mb 的远端 11q 不在其中，且仍只存在于一条 homolog 上，因此该 junction 与远端丢失不矛盾。它与 fold-back 的关系尚未确定；第 7.4 节列出三种可能的解读。 |",
])
table_notes = table_notes.replace(old_s05, new_s05).replace("3. **第 1–2 行。**", "3. **第 1–4 行。**")
assert "第 1–4 行" in table_notes

new_s2 = """## 2. 核心结论（Core Conclusions）

本节每条给出一句话结论、对客户意味着什么、证据在哪里。第 3 节按主题列出测量结果；第 8.1 节说明每条结论的依据；第 8.2 节按样本组更详细地解释每项结果；第 8.3 节说明数据不能说明什么；第 8.4 节为每个未决问题配一个实验。

1. **S04（+Dox 克隆 3）：断端已解析 —— chr11 的 114,605,216–218 连接到端粒重复，其远端一条 chr11q homolog 一直丢失到 11q 末端。** 意味着什么：该克隆带有 11q 末端缺失，其断端被端粒重复封顶。可用 junction-specific 的 PCR 引物对（第 7.11 节）在 S04 DNA 上验证，以 S01–S03 作阴性对照。证据：第 7.4 节，`tables/junction_table_target.tsv`（J0922 行），`figures/junction_views/Fig4_J0922.png`，`tables/key_junction_read_evidence.tsv`，`locus_BAM/S04_DLD1_c3_dox.locus_rescue.bam`。
2. **S05（+Dox 克隆 4）：丢失的是同一条 chr11q homolog，但断端不同且更复杂 —— 断崖处一个 376 bp 的 fold-back inversion，以及断崖远端 1.2 kb 处的一个 chr21 junction。** 意味着什么：在 fold-back 处，断端折回到自身（不是插入）；chr21 junction 只带走约 1 kb 的 chr11 序列，与远端丢失不矛盾；两个事件之间的关系尚未确定（第 7.4 节）。chr21 junction 可用 J0786 引物对验证（第 7.11 节）；fold-back 需要另行设计（第 8.4 节）。证据：第 7.4 与 7.7 节，J0359 与 J0786 行，`figures/junction_views/Fig4_J0359.png` 与 `Fig4_J0786.png`，`supporting_analyses/chr11_BAF/retained_homolog_S04_vs_S05.tsv`，`tables/caller_concordance_by_sample.tsv`。
3. **S01–S03 仍是有效基线：三者都不带三个 +Dox junction 中的任何一个。** 意味着什么：可在 junction PCR 中作阴性对照。S04 属于一个没有任何对照克隆代表的亚系（sub-lineage），因此 S04 中标为 +Dox-specific 的事件无法与该亚系在 doxycycline 处理之前已有的变异区分（第 8.3 节）；S05 不受此影响。证据：第 7.4 节表格（reference fragment 计数），`tables/junction_by_sample_long.tsv`。
4. **GM12878 克隆：位点上没有超过一个 fragment 的克隆特异重排。** 意味着什么：S07–S10 在该位点没有需要跟进的事件；目标区近端 13 kb 处、五个 GM12878 样本都有的那对 junction 是 GM12878 细胞系本身的特征，不属于任何一个克隆。证据：第 7.5 与 7.8 节，`tables/junction_table_target.tsv` 的 `GM12878_clone_specific` 列，`tables/junction_table_window_context.tsv` 的 J1365 与 J0533 行。
5. **该位点的大多数 junction 候选是 22 kb repeat 内的单 read，不应逐个跟进。** 意味着什么：527 个候选中，read rescue 后有 4 个具有两个或以上独立 fragment，其中三个就是第 1–2 条的 +Dox junction；后续工作应集中在这些上。证据：第 7.1 与 7.5 节，图 5，`tables/junction_table_target.tsv` 的 `classification` 与 `claim_level` 列。
6. **全基因组：13 个候选结构变异只属于其中一个 +Dox 克隆（S04 8 个，S05 5 个）；没有一个为两个克隆共有，也没有一个位于 chr11。** 意味着什么：没有发现两个 +Dox 克隆共有的第二处重排。两个最大的事件中，S04 的 327 kb chr8 duplication 经 read depth 与等位基因平衡确认；S05 的 chr4 事件是真实的 junction，但不是缺失；其余 11 个是候选，可在全基因组 BAM 文件中查看并做 PCR（第 8.4 节）。证据：第 7.9 节，`genome_wide_DLD1_comparison/dld1_dox_specific_events.tsv`，`supporting_analyses/SV_verification/`。

---

"""
new_p2 = ("### 8.2 实用解读（Practical reading）\n\n"
          "本节按样本组解释每项结果是什么、意味着什么，比第 2 节更详细；表后的注释说明 S04 与 S05 各行的解读限度。\n\n"
          + table_notes + "\n\n")
t = t[:s2] + new_s2 + t[s3:p2] + new_p2 + t[p3:]

# 3. Methods: UCSC table names
old = "Repeat 注释：UCSC hg38 RepeatMasker、simpleRepeat 与 genomicSuperDups 表。"
assert t.count(old) == 1
t = t.replace(old, "Repeat 注释：UCSC hg38 RepeatMasker、simpleRepeat（由 Tandem Repeats Finder 检出的串联重复）与 genomicSuperDups（片段重复，segmental duplications）表。")

# 4. Section 7.4: three readings of the chr21 junction, after the S05 paragraph
anchor = "因为 S04 已丢失其中一条拷贝的这一段。"
assert t.count(anchor) == 1
t = t.replace(anchor, anchor + "\n\n"
    "chr21 junction 与 fold-back、与远端丢失之间的关系尚未确定。有三种解读仍然成立，现有数据无法区分：（i）fold-back 复制出的那一段再次断裂并连接到 chr21，即在同一条 homolog 上发生第二轮断裂与重连；（ii）该 junction 只存在于克隆中的一部分细胞：其等位基因频率为 0.27（3 个 junction fragment 与 8 个 reference fragment；95% CI 0.10–0.57），低于每个细胞都在一条 chr11 homolog 上带有该 junction 时预期的 0.5，但只有三个 junction fragment，无法与 0.5 区分；（iii）它是与 11q 丢失无关的独立事件。未检查 chr21 的拷贝数。用 J0786 引物对（第 7.11 节）做跨 junction 的 PCR，可检验该 junction 是否真实；它位于哪条 homolog、是否每个细胞都带，需要进一步的检测，例如单细胞 PCR。")

open(CN, "w").write(t)
print("CN updated")
