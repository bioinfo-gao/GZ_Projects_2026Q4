# P42 内部分析总结（2026-10-23）

## 当前结论

- 全 18 个 WGS 样本、59 个 local assembly loci、90 条 contig、GRIDSS/SvABA external validation 与 fixed-candidate raw re-genotyping 均已完成。
- 282 个 external-caller BND 中，94 个可解析为双端几何；4 个仅在 DLD1 `+dox` cases 有 raw dual-endpoint QNAME evidence，但没有一个达到预注册的 Dox/clone-specific high-confidence gate。
- 四个 case-only candidates 的 S01/S02/S03 standard local depth 为 14.7–23.1×，低于 `3 fragments / 5%` 所需的 60×；C2 的 S05 三个 raw QNAME 均不属于 standard-eligible subset。所有四个也缺少 sequence-level contig junction。
- 因此对外允许的结论是：存在四个明确保存的 case-only raw-read candidates，但短读长 WGS 本轮**未建立**高置信 Dox/clone-specific breakpoint；不得写成“controls 中不存在”或“FISH 已被 WGS 证实”。

## 关键审计数

| 项目 | 数值 | 解释 |
| :--- | :---: | :--- |
| 输入 WGS samples | 18 | S01–S18，样本表为 `scripts/0_sample_table.csv` |
| Local assembly loci | 59 / 59 completed | 90 条 contig；32 条有 alternative loci |
| Normalized GRIDSS events | 278 | 18-sample caller run；caller fields 不是 independent fragments |
| SvABA germline / somatic events | 4 / 0 | corrected five-sample DLD1 invocation |
| Fixed caller BND universe | 282 | 包含所有 stage-7 normalized BND events |
| Two-endpoint geometries | 94 | 188 single-anchor/unparsed 保留，不能视作阴性 |
| Case-only raw candidates | 4 | C1–C4；均未达到 final gate |
| Raw evidence rows, S01–S05 | 51,508 | 19,557 anchor-context matches；31,951 未链接到 fixed caller geometry，仍留在 discovery universe |

## 当前有效产物

- Current delivery figure suite: `figures/20261023_fullcohort_repeat_resolve_delivery_v7/`（12 PNG + 12 PDF；v1–v6 为已保留的试运行/非最终版本）。
- Fixed-candidate source of truth: `analysis/20261023_fullcohort_repeat_resolve/09_fixed_candidate_regenotyping_v1/`。
- Prepared, not sent, client package: `custom_research_report_20261023/`。

## 仍需用户过目的边界

1. 客户交付报告将 FISH 写为 qualitative observation，因为分析输入没有 raw FISH image/table；报告不将其量化也不把它与 WGS 合并计数。
2. 建议的后续 confirmation 是针对 C1–C4 的 breakpoint-spanning targeted PCR；若后续 WGS，建议提高 local depth。没有提出 long-read platform。
3. 本轮交付尚未发送；用户过目并明确允许后才可对外发送。
