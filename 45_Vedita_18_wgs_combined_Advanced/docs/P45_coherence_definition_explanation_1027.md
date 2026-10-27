# P45「一致性（coherent with the copy-number step）」：精确含义与报告改写建议

- 日期：2026-10-27
- 状态：内部文档，待审阅。审阅完毕后，与其他已批准的改动一起进入下一版报告源文件和交付文件夹；本文档不改动任何报告。
- 涉及位置：报告源文件 `docs/_client_report_source.md` 第 51、124–127、207–208、367、397 行；交付表 `custom_research_report_20261026/tables/structural_coherence_tests.tsv`。

## 0. 一句话

**"一致性"是在问：这个 junction 是否正好位于拷贝数下降的位置，而且下降的那一侧确实丢了一条拷贝。** 两条都满足，junction 读段和拷贝数、BAF 就讲的是同一件事："chr11 在这里断开，断点以远少了一条"。

## 1. 精确定义：两项检验，同一样本内同时满足

只对在某样本中有 ≥ 2 个 fragment 支持的 junction 做检验，所有数据都取自该样本。

| 检验 | 精确条件 | 回答的问题 | 列名（交付表） |
| :--- | :---: | :---: | :---: |
| 检验 1：位置吻合 | junction 的某个 chr11 断点区间（含 microhomology 可滑动范围）与该样本深度断崖边界的 95% CI 有重叠 | junction 是不是在拷贝数发生变化的地方 | `criterion1_in_step_ci` |
| 检验 2：确有丢失 | 在 S01–S03 中都杂合的位点，断崖远端（114.7–135.0 Mb）有 ≥ 90% 在该样本中变成纯合（BAF < 0.1 或 > 0.9） | 断崖以远是否真的只剩一条拷贝 | `criterion2_LOH_ge_0.90` |
| 一致性 | 检验 1 和检验 2 同时为 1 | — | `coherent` |

## 2. 为什么两项都要

**只看检验 1 不够。** 没有真实断崖的样本，断崖边界的 95% CI 会宽到覆盖整个区域，任何 junction 都会"落在 CI 内"。交付表里正是这样：

| junction | 样本 | 深度断崖 95% CI | 检验 1 | 远端 LOH 比例 | 检验 2 | 一致性 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| J0922 | S04 | 114,605,000–114,605,599（600 bp） | 1 | 0.978 | 1 | **1** |
| J0359 | S05 | 114,605,600–114,606,199（600 bp） | 1 | 0.98 | 1 | **1** |
| J1635 | S01 | 114,592,800–114,631,199（38 kb） | 1 | 0.0 | 0 | 0 |
| J1635 | S02 | 114,590,600–114,639,199（49 kb） | 1 | 0.0 | 0 | 0 |
| J1635 | S03 | 114,592,200–114,636,599（44 kb） | 1 | 0.0 | 0 | 0 |
| J1036 | S03 | 114,592,200–114,636,599（44 kb） | 1 | 0.0 | 0 | 0 |

J1635、J1036 在基线样本里"通过"了检验 1，只是因为基线没有断崖、CI 覆盖了整个窗口；检验 2 把它们排除了。所以：**检验 2 证明断崖真实存在，检验 1 证明 junction 在断崖上。** 交付表共 31 行，一致性为 1 的只有 J0922（S04）和 J0359（S05）两行。

**J0786（S05 的 chr21 junction）为什么不算**：它位于 114,607,190–192，在 S05 断崖 CI（114,605,600–114,606,199）远端约 1 kb，检验 1 不满足，所以停在 "limited junction candidate"。

## 3. 一致性能说明什么、不能说明什么

| 能说明 | 不能说明 |
| :--- | :---: |
| junction 与拷贝数丢失在位置上吻合，是丢失的合理断端 | 它就是造成丢失的那个断端（同一位置可能还有未检出的其他结构） |
| 三种信号（junction 读段、depth、BAF）互相印证 | 独立验证：三者都来自同一次测序 |
| 可以把 partner 端不明的 junction（J0922，接端粒重复）升一级 | partner 端在哪里：一致性只检验 chr11 端 |

所以最高等级叫 "limited … coherent with the copy-number step"，而不是 confirmed；confirmed 需要 PCR 或 FISH。

## 4. 那句"除外"的完整意思（第 124 行）

- **一般规则**：junction 只要有一端不能唯一定位，就是 "Repeat-ambiguous junction"，不作进一步解读。
- **例外**：不明的只是 partner 端（chr11 端唯一确定），并且通过上面两项检验，就升到 "coherent with the copy-number step"。
- **理由**：两项检验只用 chr11 端的位置；partner 端接到哪里不影响"chr11 在这里断开、以远丢了一条"这个判断。
- **适用实例**：只有 J0922（S04）。chr11 端在 114,605,216–218 唯一确定；partner 是 (TTAGGG)n，chrX、chr2、chr4、chr10 都能同样好地比对。
- **不适用**：chr11 端本身不明（例如断点落在 22 kb tandem array 里、能在重复拷贝间滑动），连检验 1 都无法做。

## 5. 现有措辞的问题

| 位置 | 现在的写法 | 问题 |
| :--- | :---: | :---: |
| 第 124 行 | "except when only the partner end is ambiguous and the event is coherent (next rows)" | "coherent" 未定义就使用；"(next rows)" 只指方向，读者要自己去下面几行拼出含义 |
| 第 125–126 行 | "not coherent" | 同上 |
| 第 127 行 | 两项条件写在一个长句里 | 没说明为什么需要两项；没说明基线样本的 CI 会覆盖整个窗口 |
| 第 207–208 行 | "coherent with step" | 表格里的简写，正文没有对应解释 |
| 第 367 行 | "the copy-number and BAF coherence of Section 5" | 又换了一种说法 |
| 第 397 行 | 交付表说明只列内容 | 没说明 `coherent` 列怎么由两列得出 |

## 6. 建议的报告改法（英文，待你审阅）

**(a) 在 §5 Claim level 表之前加一段定义（新增）：**

> **Coherence with the copy-number step.** A junction is called *coherent* in a sample when two tests, both computed in that sample, are met: (1) *position*: one of its chr11 breakpoint intervals overlaps the 95% confidence interval of that sample's depth step; and (2) *loss*: at least 90% of the sites heterozygous in all of S01–S03 are homozygous in that sample distal of the step (114.7–135.0 Mb), i.e. one copy of distal 11q is really lost. Test (2) is needed because a sample without a step has a confidence interval spanning the whole window, so test (1) alone would be met by any junction. Coherence uses only the chr11 end of the junction; it says that the junction lies where one copy is lost, not that it is the only structure there, and it is not independent confirmation, because the junction reads, the depth and the allele fractions all come from the same sequencing run. The tests were applied to every junction with two or more fragments in a sample (`tables/structural_coherence_tests.tsv`).

**(b) 第 124 行改为：**

> Ambiguous (repeat). Exception: if only the partner end is ambiguous and the chr11 end is unique, the junction can still be coherent (definition above) and is then classed in the last-but-one row

**(c) 第 125–126 行**："not coherent" 保持，前面已有定义。

**(d) 第 127 行改为：**

> Limited, Well supported, or Ambiguous only at the partner end; and coherent with the copy-number step in a sample with ≥ 2 fragments (both tests of the definition above)

**(e) 第 367 行**：把 "the copy-number and BAF coherence of Section 5" 改为 "the coherence tests of Section 5 (depth step and allele fractions)"，与定义用同一术语。

**(f) 第 397 行**，交付表说明末尾加：

> ; column `coherent` is 1 when `criterion1_in_step_ci` and `criterion2_LOH_ge_0.90` are both 1

**(g) 第 207–208 行表格**：保持 "coherent with step"，正文已有定义。

以上英文措辞都是建议，需要你审阅和修改。

## 7. 复核

```
awk -F'\t' 'NR==1 || $9==1' /home/gao/projects_2026Q4/45_Vedita_18_wgs_combined_Advanced/custom_research_report_20261026/tables/structural_coherence_tests.tsv
```

2026-10-27 实测输出：表头加 J0359（S05）、J0922（S04）两行。
