# P45 下一版交付文件夹的待办（审阅完毕后集中生成）

- 用户指示（2026-10-27）：整个报告审阅完毕后，集中生成下一版交付文件夹；在此之前不生成新文件夹。
- 生成命令：`bash scripts/20261025/12_build_delivery.sh <生成当天日期>`；1025、1026 文件夹不动。

| # | 事项 | 状态 | 依据 |
| :--- | :---: | :---: | :---: |
| 1 | chr4 "888 kb deletion" 改为无拷贝丢失的接头；chr8 327 kb 获得注明已核实（第 7a 条、§7.9 两行与表后说明） | 用户已同意，已写入报告源文件（提交 `96b031f`） | `docs/P45_dox_specific_SV_chr4_chr8_verification_1027.md` |
| 2 | 第 1 项引用的核查数据（`analysis/20261027_sv_verification/`）是否进交付表格 | 待用户决定 | 同上 |
| 3 | "同一条 chr11q homolog" 的明确表述（Key Finding 4、第 40、41、281、381 行） | 待审阅 | `docs/P45_same_chr11q_homolog_explanation_1027.md` |
| 4 | "一致性（coherent with the copy-number step）"的定义与改写（§5 新增定义段，第 124、127、367、397 行） | 用户已同意（2026-10-27），已写入报告源文件 | `docs/P45_coherence_definition_explanation_1027.md` |
| 5 | 下一版中文版：复制成新日期文件后再改，1025、1026 中文版不动；同步第 1、3、4 项 | 构建时做 | — |
| 6 | 图 1 图注下层 panel 加注释：底部 y = 1 的密集 Weak 点是各样本共有的本底，可作对照 | 待审阅（措辞见下） | `custom_research_report_20261026/tables/junction_table_target.tsv`；`scripts/20261025/11_figures.py` 第 143–159 行 |

## 第 6 项：图 1 下层 panel 注释

**数据（2026-10-27 由交付表 `junction_table_target.tsv` 计算）**

- 下层 panel 共画 258 个 junction（任一样本 ≥ 1 个 rescue 后 fragment）；其中 254 个是 Weak（1 个 fragment），即 y = 1 那一排。
- 254 个 Weak 中 247 个只出现在一个样本里，6 个出现在 2 个样本里，1 个出现在 4 个样本里：不重复出现。
- 每个样本的 ≥ 1 fragment junction 数：S01–S03 为 27 / 21 / 26（其中 Weak 27 / 21 / 25）；S06–S10 为 29 / 20 / 36 / 40 / 29（全部 Weak）；S04、S05 为 22 / 17（其中 Weak 21 / 15）。没有深度断崖的样本密度相同，+Dox 克隆并不更多。
- 带 S04-only / S05-only 标记的点共 39 个（标记规则见 `11_figures.py` 第 152 行：`case_specific` 与 `case_specific_underpowered` 都画标记）；其中 36 个是 Weak（S04 21 个、S05 15 个），只有 3 个有 ≥ 2 个 fragment：J0922（Ambiguous_repeat，S04）、J0359 和 J0786（Limited，S05）。

**解读**：在没有重排的样本里也有同样密度的单 fragment junction，说明这一排是 tandem repeat 内 read 错误定位产生的本底，是一个内置的阴性对照；单个 fragment 的 "+Dox-specific" 标记不高于这个本底，不作解读。高出本底的只有 3 个 ≥ 2 fragment 的 junction。

**建议加在英文图注 "marker shape = +Dox-specific events." 之后（待审阅）**：

> Most points (254 of 258) are Weak junctions with a single fragment, and 247 of these occur in only one sample. They appear at a similar density in every sample, including those without a copy-number step: 21–27 per DLD1 baseline sample (S01–S03), 20–40 per GM12878 sample (S06–S10) and 15–21 in the +Dox clones (S04, S05). This row at the bottom is therefore an internal control: the level of single-fragment, junction-like reads that mis-placement within the repeat produces without any rearrangement. Of the 39 points marked as S04-only or S05-only, 36 are such single-fragment points and are not interpreted; only J0922, J0359 and J0786, with two or more fragments, stand above this background.
