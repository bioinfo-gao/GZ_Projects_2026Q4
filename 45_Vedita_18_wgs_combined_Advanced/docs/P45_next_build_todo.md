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
