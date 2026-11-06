# evals/ —— 哪个是当前有效版

本目录**只放提交给 LatchBio 的交付物**，共两个 JSON，一题一个。

| 文件 | 状态 | 说明 |
| :--- | :---: | :--- |
| `bbb_label_implication_v1.json` | ✅ **当前有效，交付** | Task 1，NumericTolerance，7 字段 |
| `glioma_margin_rate_limit_v4.json` | ✅ **当前有效，交付** | Task 2，MultipleChoice，5 字段 |

## Task 2 的版本沿革（四代，均因泄题退役）

| 版本 | 退役原因 | 泄题层级 | 归档 |
| :--- | :--- | :---: | :--- |
| v1 | 前提写 `growth-arrested`、字段名含 `arrest`，与选项 `quiescence` 近义 | 词形 | `superseded_for_comparison/evals/v1_20261105_lexical_leak_growth_arrested/` |
| v2 | "选择因三药共有属性而成立的因素" ⇒ 给排除法立路标 | 解法 | `superseded_for_comparison/evals/v2_20261105_elimination_solvable/` |
| v3 | 四条"临床事实"陈述了 GBM 一般生物学，第 1 条恰是让 BBB 成为诱骗项的那条知识 | 知识 | `superseded_for_comparison/evals/v3_20261105_premises_stated_the_biology/` |
| **v4** | **当前有效** | — | — |

**三代共有同一个错误：把本该由作答者完成的那一步替他做了。**
v4 因此新增了显式规格原则（见其 `notes.specification_principle`）：

> **患者层面、因病例而异、无法先验推知的 ⇒ 写明。**
> **疾病层面、对所有 GBM 都成立的 ⇒ 一律不写，那正是考题。**

废弃的 JSON 按用户指示（2026-11-05「过时的文件转入 tmp」）移至 `tmp/`，
逐字节副本留在上表的归档目录中，因此移动不丢失唯一副本。
⇒ **本目录内不存在非交付文件。**
