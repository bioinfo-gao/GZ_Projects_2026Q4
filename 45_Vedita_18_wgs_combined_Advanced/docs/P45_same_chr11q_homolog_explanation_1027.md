# P45「两个 +Dox 克隆丢失的是同一条 chr11q homolog」：含义、证据与报告措辞

- 日期：2026-10-27
- 状态：内部文档，待审阅。审阅完毕后，与其他已批准的改动一起进入下一版报告源文件和交付文件夹；本文档本身不改动任何报告。
- 涉及的报告位置：报告源文件 `docs/_client_report_source.md` 第 40、41、69、281、348、349、358、381 行（见 §5）。

## 0. 结论

**两个克隆丢的是同一条。** 在 chr11 约 114.6 Mb 以远（到 q 臂末端），S04 和 S05 丢失的是两条 chr11 中的**同一条**，保留下来的也是同一条。这不是"两个克隆各丢了一条，但不知道是哪条"。

**确定不了的是这条来自父亲还是母亲。** 没有 DLD1 供体父母的基因型，所以这一条只能称为"同一条 homolog"，不能标为父源或母源。

## 1. 术语

| 用语 | 本文含义 |
| :--- | :---: |
| homolog（同源染色体） | 二倍体细胞中两条 chr11 之一，一条来自父亲，一条来自母亲 |
| 杂合位点 | 两条 homolog 在该位置带不同碱基（等位基因 A / B）的位点；本文只用在 S01–S03 中都杂合的位点 |
| BAF（B-allele frequency，B 等位基因频率） | 该位点上带 alt 等位基因的 reads 比例；两条都在时约 0.5，只剩一条时接近 0 或 1 |
| LOH（loss of heterozygosity，杂合性丢失） | 原本杂合的位点在克隆中变成纯合（报告的判据：BAF < 0.1 或 > 0.9） |
| "同一条 homolog" | 两个克隆剩下的那条 chr11q 来自同一亲本拷贝。由于两个克隆都来自同一个 DLD1 群体（S01），它们的两条 chr11 分别是同一对亲本 homolog 的后代拷贝，所以可以逐位点比较 |

## 2. 证据与推理

### 2.1 每个克隆都只剩一条

| 样本 | 远端（114.7–135.0 Mb）可评估的杂合位点 | 其中变为 LOH 的比例 | 远端 depth |
| :--- | :---: | :---: | :---: |
| S01–S03 | 16,650 | 0.000 | 基线 |
| S04 | 10,464 | **0.978** | 约 0.50× |
| S05 | 12,830 | **0.980** | 约 0.51× |

文件：`custom_research_report_20261026/supporting_analyses/chr11_BAF/baf_segments.tsv`。depth 数据见报告第 7.4–7.6 节。

**depth 减半，同时约 98% 的杂合位点变成纯合，说明这一段只剩一条拷贝**，而不是复制出两份同一条（后者 depth 不会减半）。

### 2.2 两个克隆剩下的是同一条

| 区段 | 两个克隆都 LOH 的位点 | 两个克隆保留同一等位基因 | 比例 |
| :--- | :---: | :---: | :---: |
| 114.7–135.0 Mb | 8,004 | 7,986 | **99.78%** |

文件：`custom_research_report_20261026/supporting_analyses/chr11_BAF/retained_homolog_S04_vs_S05.tsv`。

两种假设下的预期一致率：

| 假设 | 预期一致率 | 实测 |
| :--- | :---: | :---: |
| 两个克隆丢的是**同一条** homolog | 约 100%（只受基因分型误差影响） | **99.78%** |
| 两个克隆丢的是**不同**的 homolog | 约 0%（一个处处剩 A，另一个处处剩 B） | — |

- 一致率在两种假设下分别接近 100% 和 0%，不存在居中的结果，所以这项检验没有歧义。
- 不一致的 18 个位点占 0.22%。它们不影响结论（见 §2.4），但其成因未查明：2026-10-27 补查发现：这 18 个两克隆之间不一致的位点，全部是 S04 保留 ref（BAF 0.00）、S05 保留 alt（BAF 1.00，其中 1 个为 0.91），反向组合（S04 保留 alt、S05 保留 ref）为 0 个；随机分型误差应两种组合大致各半，所以不像随机误差。此前写的「与分型误差同量级、不需另找解释」已收回。

### 2.3 两次丢失是独立发生的

两个克隆的断端不同：

| 克隆 | 断端 | 报告位置 |
| :--- | :---: | :---: |
| S04 | chr11 序列止于 114,605,216–114,605,218，直接接端粒重复 | 第 40、348 行 |
| S05 | 114,605,973 / 114,606,349 处的 fold-back inversion，以及 1.2 kb 以远的 chr21 接头 | 第 41、349 行 |

断端结构不同，说明不是同一个祖先细胞里的一次事件被两个克隆共同继承，而是**两次独立断裂，打在了同一条 homolog 上**。

### 2.4 一致性沿整个远端逐段成立（2026-10-27 补算）

如果两个克隆只在部分区段丢的是不同的 homolog（例如其中一个克隆在远端某处发生过有丝分裂重组），平均 99.78% 会掩盖这一段：那一段的一致率会接近 0%。按 2 Mb 分段重算（数据：`baf_sites.tsv.gz`，判据同 §2.2）：

| 区段（Mb） | 两克隆都 LOH 的位点 | 保留同一等位基因 | 一致率 |
| :--- | :---: | :---: | :---: |
| 114.7–116.7 | 819 | 819 | 1.0000 |
| 116.7–118.7 | 602 | 602 | 1.0000 |
| 118.7–120.7 | 699 | 699 | 1.0000 |
| 120.7–122.7 | 601 | 599 | 0.9967 |
| 122.7–124.7 | 800 | 799 | 0.9988 |
| 124.7–126.7 | 894 | 892 | 0.9978 |
| 126.7–128.7 | 674 | 667 | 0.9896 |
| 128.7–130.7 | 665 | 664 | 0.9985 |
| 130.7–132.7 | 1,039 | 1,035 | 0.9962 |
| 132.7–134.7 | 988 | 987 | 0.9990 |
| 134.7–135.0 | 223 | 223 | 1.0000 |

- 11 段全部 ≥ 98.96%，没有任何一段接近 0%：从断崖到 11q 末端，两个克隆剩下的都是同一条。
- 18 个不一致位点散布在 120.93–134.43 Mb，最密的 126.7–128.7 Mb 段有 7 个（674 个位点中），不构成一段连续的相反等位基因。
- 这 18 个不一致位点全部是同一种组合（S04 保留 ref、S05 保留 alt；反向组合 0 个），成因未查（可能是位点级比对偏差，也可能是某一克隆的点突变，未核实）；无论哪种，都不改变「同一条 homolog」的结论。

## 3. 能说什么、不能说什么

| 说法 | 是否可以 | 依据 |
| :--- | :---: | :---: |
| 两个克隆在 114.6 Mb 以远丢失的是同一条 chr11q homolog | 可以 | §2.1 + §2.2 |
| 两个克隆保留的是同一条 homolog | 可以（与上一条等价） | §2.2 |
| 这是两次独立事件 | 可以 | §2.3 |
| 丢失的是父源 / 母源 chr11 | **不可以** | 没有父母基因型，也没有跨亲本的定相信息 |
| 该 homolog 对 +Dox 诱导的断裂更敏感（存在偏好） | **不可以** | 两个独立克隆碰巧丢同一条的概率本来就是 1/2，n = 2 不足以检验偏好 |
| 近端（114.6 Mb 以内）也是同一条 homolog 发生变化 | **未检验** | 报告只比较了远端；S05 近端的渐变失衡（BAF 偏离 0.159）发生在哪一条上，没有做逐位点比较 |

**关于偏好的说明**：要检验是否偏向这一条，需要更多独立的 +Dox 克隆。例如 6 个独立克隆全部打在同一条上，在无偏好假设下的概率为 2 × (1/2)^6 ≈ 3%。这正是报告 Next Steps 第 381 行建议的内容。

## 4. 现有报告措辞的问题

| 位置 | 现在的写法 | 问题 |
| :--- | :---: | :---: |
| 第 69 行 Key Finding 4 | "Both +Dox clones lost the same chr11q homolog." | 没写区段；没说明依据；读者可能理解为"都丢了一条 chr11q"；没说明无法区分父源/母源；没提醒偶然概率为 1/2 |
| 第 40 行（S04） | "the distal chr11q copy that is lost is one complete homolog" | "one complete homolog" 容易读成"整条 chr11"，实际只是 114.6 Mb 以远的那一段 |
| 第 41 行（S05） | "The same distal homolog is lost as in S04" | 含义正确，但没有给出依据 |
| 第 281 行（第 7.7 节） | "the same one in both (8,004 sites homozygous in both clones; the same allele retained at 99.78%)" | 有数字，但没有说明这个数字为什么能区分"同一条"和"不同条" |
| 第 348、349 行（结论表） | "Terminal loss of one chr11q homolog" / "Loss of the same homolog" | 可保留；只需与 Key Finding 4 的用语一致 |
| 第 358 行（Conclusion 2） | "on the same homolog" | 同上 |

## 5. 建议的报告措辞（英文，待你审阅）

| 位置 | 建议改为 |
| :--- | :---: |
| 第 69 行 Key Finding 4 | **Both +Dox clones lost the same one of the two chr11 homologs distal of 114.6 Mb.** Evidence column: "At 7,986 of 8,004 sites (99.78%) that became homozygous in both clones, the two clones retained the same allele; had they lost different homologs, they would retain opposite alleles at nearly every site. The two losses have different broken ends and are therefore independent events. Whether this homolog is maternal or paternal cannot be determined without parental genotypes; with two independent clones, the same homolog is expected by chance in half of cases. Section 7.7." |
| 第 40 行（S04） | 把 "the distal chr11q copy that is lost is one complete homolog" 改为 "distal of the break, one of the two chr11 homologs is lost to the end of 11q" |
| 第 41 行（S05） | 把 "The same distal homolog is lost as in S04" 改为 "Distal of the step, the same homolog is lost as in S04 (Section 7.7)" |
| 第 281 行（第 7.7 节） | 在 "...the same allele retained at 99.78%)" 后加："If the two clones had lost different homologs, they would retain opposite alleles at almost every one of these sites, so agreement near 100% identifies the same homolog; the agreement holds along the whole segment (≥ 98.96% in each of eleven 2 Mb windows from the step to the end of 11q), so no part of the distal arm was lost from different homologs in the two clones; the 18 discordant sites (0.22%) are scattered single sites, not a contiguous segment. Which parent this homolog came from is not known. The comparison covers the distal segment only; on which homolog the graded proximal imbalance of S05 lies was not tested." |
| 第 348、349、358 行 | 不改 |
| 第 381 行（Next Steps） | 在问题后加："With two clones, the same homolog is expected by chance in half of cases; a preference would require more independent clones (for example, all of six clones on the same homolog has a probability of about 3% by chance)." |

以上英文措辞都是建议，需要你审阅和修改。

## 6. 复核方法

- 两个克隆都 LOH 的位点及一致率：`custom_research_report_20261026/supporting_analyses/chr11_BAF/retained_homolog_S04_vs_S05.tsv`（1 行），由 `scripts/20261025/13_self_audit.py` 第 158 行起与报告文字核对。
- 各区段的 LOH 比例：同一目录下的 `baf_segments.tsv`。
- 逐位点 BAF：同一目录下的 `baf_sites.tsv.gz`，175,122 行数据；第 1 列为细胞系，DLD1 行的列依次为 `pos` 和 S01–S05 的 BAF。下面的命令可重算一致率，2026-10-27 实测输出 `8004 7986 0.997751`，与交付表格一致：

```
zcat baf_sites.tsv.gz | awk -F'\t' 'NR>1 && $1=="DLD1-136" && $2>114700000 && $2<=135000000 && $6!="" && $7!="" && $6!="NA" && $7!="NA" { a=($6<0.1)?"R":($6>0.9)?"A":""; b=($7<0.1)?"R":($7>0.9)?"A":""; if(a!="" && b!=""){n++; if(a==b) s++} } END{print n, s, s/n}'
```
