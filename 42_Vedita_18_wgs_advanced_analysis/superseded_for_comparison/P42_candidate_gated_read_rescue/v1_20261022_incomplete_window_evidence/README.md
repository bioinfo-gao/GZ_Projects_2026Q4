# v1 — candidate-gated read-rescue（保留作比较，不作当前主线结论）

- **归档日期**：2026-10-23
- **原始运行日期**：2026-10-22
- **归档原因**：该轮从 P41 已保留的 179 个 candidate 及其 source read-name 开始，并未建立目标窗口的全量证据宇宙。因此它不能回答 repeat-rich chr11 区域中真实定位、遗漏 junction、secondary/multi-mapping 解释、unmapped mate 或 clone-specific occurrence 的问题。
- **保留价值**：它仍是对 179 个 P41 预选候选所做的 source-sequence re-alignment 记录，可与新流程逐 candidate 比较；不得删除、覆盖或伪称为全窗口 discovery。

## 归档内容

| 路径 | 内容 |
| :--- | :---: |
| `analysis/20261022/` | 旧 candidate table、source-sequence extraction、re-alignment 与汇总结果 |
| `scripts/20261022/` | 产生该版本的完整脚本快照 |
| `scripts/0_sample_table_10_samples.csv` | 旧 10-sample SSOT 快照 |
| `delivery/custom_research_report_supplementary_20261022/` | 当时准备但未发送的 supplementary delivery |
| `docs/` | 当时的计划、方法/结论和问题审计文档快照 |

## 与新主线的关系

2026-10-23 起的新主线使用 **18 个样本**、全窗口 read inventory、逐分子而非逐 record 计数、local assembly、depth/CNV breakpoint 定位、Manta/DELLY caller concordance 与 DLD1 clone-specific 比较。新输出只写入新的带日期目录，绝不原地改写本版本。
