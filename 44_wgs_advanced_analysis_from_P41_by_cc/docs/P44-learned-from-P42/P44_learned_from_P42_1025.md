# P44 向 P42 学习：交付对比评估与改进计划

- **创建：** 2026-10-25（文件名日期 = 创建日，永不改名；修订记在文末更新记录）
- **性质：** 内部文档，不进客户交付。
- **评估对象：** P42（Codex Agent）交付目录 `/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis/custom_research_report_20261024/`，报告 `P42_Vedita_18_WGS_1024.md`。
- **对照对象：** P44 交付目录 `/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc/custom_research_report_supplementary_20261025/`（生成日 2026-10-25，**尚未经用户审阅，未发送**）。
- **相关文档：** 研究计划 [`../P44_research_plan_and_review_1024.md`](../P44_research_plan_and_review_1024.md)；运行记录 [`../P44_key_operations_and_parameters_log_1025.md`](../P44_key_operations_and_parameters_log_1025.md)（其 §4 已记 P42 C1–C4 与 P44 的坐标对照，本文不重复）。

---

## 0. 结论先行

1. **P42 在交付工程上强于 P44**：header 去本机路径并逐文件审计、客户问题 → 答案 → 证据文件的机读对照表、逐文件数据清单、外部输入缺口表、关键候选的 in-silico 引物与全基因组特异性筛查、关键事件逐 read 证据账本。这些 P44 都没有或只有文字版。
2. **P44 在科学产出与客户可读性上强于 P42**：P42 把分析收窄到 C1–C4 四个候选，漏掉了本项目最关键的两个发现（S04 的 chr11→(TTAGGG)n 端粒 junction、S05 的 376 bp fold-back inversion），客户的"rescue 前后逐事件计数"也只回答了 4 个候选。报告大量内部审计术语，§2 没有给出"对客户意味着什么"。
3. **P44 有一个必须修的缺陷**：交付的 BAM/VCF header 带我方本机路径（实测见 §3.1）。另有 5 项值得借鉴的交付做法（§4）。
4. **不借鉴**：P42 的候选封闭式分析框架、审计口吻的客户正文、客户文字里的内部项目编号、把 846 GB BAM 放进报告包（§5）。

---

## 1. 评估方法（读了什么、怎么判断）

| 读取内容 | 路径（相对 P42 交付目录） | 用途 |
| :--- | :---: | :---: |
| 报告全文（149 行） | `P42_Vedita_18_WGS_1024.md` | 结构、结论、措辞、客户可读性 |
| 目录结构与体积 | `find -L . -type f`；`du -shL .` = 846 G | 交付组织方式 |
| `results/` 20 个 TSV 的表头与行数 | `results/*.tsv` | 机读产物覆盖面 |
| `methods/` 6 个文件 | `methods/*` | 方法透明度、header 审计 |
| 图件 manifest 与 3 张代表图 | `figures/figure_manifest.json`；Fig10、Fig13、Fig14 | 图件质量 |
| BAM 是复制还是链接 | `stat -c '%h %s %n' genome_wide_data/bam/*.bam` → link count = 3 | 磁盘占用、交付方式 |

P44 一侧的对照检查（2026-10-25 实测）：

```bash
D=/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc/custom_research_report_supplementary_20261025
E=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin
$E/samtools view -H $D/locus_BAM/S04_DLD1_c3_dox.locus_rescue.bam | grep -c '/home/\|/mnt/\|/Work_bio'            # 6
$E/bcftools view -h $D/genome_wide_SV_VCF/manta/S04_DLD1_c3_dox.manta.diploidSV.PASS.vcf.gz | grep -c '/home/\|/mnt/\|/Work_bio'   # 4
$E/bcftools view -h $D/supporting_analyses/chr11_SNV/chr11_S01_S10.QUAL30.vcf.gz | grep -c '/home/\|/mnt/\|/Work_bio'           # 7
```

**判断原则**：科学上以研究计划 §2A（先定义 data universe 再定义 candidate universe、filter-impact 审计、不跨证据层级下结论）为尺；交付上以 report / delivery skill（§1–§9 结构、§2 与 §8.2 一简一详、缩写首次全称、不出现内部编号、交付必给清单）为尺。

---

## 2. P42 的优点（值得学）

| # | 优点 | 具体体现（P42 交付内路径） | 为什么有价值 |
| :--- | :--- | :---: | :--- |
| A1 | **交付文件去本机路径，并逐文件审计** | `methods/bam_header_sanitization_audit.tsv`（18 行：源/交付字节数、alignment record 数前后一致、去掉的 header 行数）；`methods/vcf_header_sanitization_audit.tsv`（90 行：variant body SHA-256 前后一致） | 我方路径（含用户名、目录结构）不外泄；同时**证明**改 header 没有动数据本体 |
| A2 | **客户问题 → 答案 → 证据文件的机读对照** | `results/customer_request_coverage.tsv`（问题、完成状态、新证据文件、边界）；`results/claim_evidence_matrix.tsv`；`results/customer_question_summary.tsv` | 客户和我们都能逐条核对"问的都答了没有"；比只在报告 §8.2 写文字更易审 |
| A3 | **逐文件数据清单** | `results/data_inventory.tsv`（109 行：类别、交付内相对路径、bytes、索引路径、availability） | 传输规划、核对完整性；与 md5 互补 |
| A4 | **外部输入缺口明示** | `results/external_input_status.tsv`：独立 PCR、敏感度研究、原始 FISH 图像"not_provided"，并写明各自对解读的影响 | 把"没做/没给"写成表，避免读者把未检验误读为阴性 |
| A5 | **"不存在的东西不说成存在"** | 报告 §9 最后一行：全基因组 SNV/indel VCF "Unavailable; not represented as present"；Manta candidateSmallIndels 不冒充 SNV/indel callset | 符合全局规则"没做过的检查不得写成带结论意味的句子" |
| A6 | **关键候选的 in-silico 引物 + 全基因组特异性筛查** | `results/assay_design_selection.tsv`、`results/primer_specificity.tsv`（46 行，每条引物全长精确命中数）、`methods/virtual_junction_templates.fasta` | §8.4 的"建议 PCR 验证"变成可直接下单的方案 |
| A7 | **关键事件逐 read 证据账本** | `results/c2_read_end_junction_resolution.tsv`：每条 read end 的两段 FLAG、坐标、CIGAR、MAPQ、query 区间、overlap 序列；`results/raw_support_alignment_audit.tsv` 带 SEQ 与 SA | 任何人可离开我们的代码独立复核 microhomology 与断点 |
| A8 | **分母口径明确区分** | `results/standard_evidence_autosomal_depth.tsv`：把"严格过滤后的分析分母深度"与"WGS 数据量 25–30×"分开定义 | 防止把过滤后深度误读为测序产出不足 |
| A9 | **DELLY BCF→VCF 转换有审计** | `methods/delly_conversion_audit.tsv`（36 行） | 格式转换可追溯 |
| A10 | **图件 manifest** | `figures/figure_manifest.json` 列 14 张图 PNG/PDF | 交付脚本可据此断言图件齐全 |

---

## 3. P42 的缺点（不要学）

### 3.1 科学层面

| # | 缺点 | 证据 | 后果 |
| :--- | :--- | :--- | :--- |
| B1 | **最终交付仍是候选封闭式**：从 GRIDSS/caller 的 282 个候选出发，只保留 4 个 case-only 候选（C1–C4）详细审 | 报告 §5："From 282 candidate observations, 94 were parseable … four case-only candidates (C1-C4) were retained" | 违反 §2A"先定义 data universe 再定义 candidate universe"；caller 没报的 junction 进不了分析 |
| B2 | **漏掉本项目最关键的两个发现** | P44 交付 `tables/junction_table_target.tsv` 中 P44J0332（S04，chr11:114,605,216–218 → (TTAGGG)n，标准计数 0 → rescue 后 4 fragment）与 P44J0115（S05，376 bp fold-back，位于 S05 深度断崖 95% CI 内）；P42 报告与 `results/` 中均无对应 | 客户最关心的"+Dox 克隆断裂末端是什么"没被回答；P42 结论变成"C1–C4 均不成立" |
| B3 | **"rescue 前后计数"只回答了 4 个候选**，且 rescue 只做了 unmapped-mate 一个来源（全 locus 回收 5 个 QNAME） | 报告 §6.5 | 客户原话要求"junction-spanning read counts per event, before and after read rescue"；XA 多重比对这一最大来源（P44 实测每样本 3,920–7,946 read pair）未纳入 rescue 计数 |
| B4 | **没有解释 chr2:32.916 Mb 吸附区的成因** | P42 报告未提 | P44 实测为短 insert 读穿 adapter 后的 poly-G 尾比对到 chr2 纯 poly-G 区，占 P41 目标区 1,181 事件中的 693 个 |
| B5 | **未发现目标区本身是 22 kb tandem repeat array**（42 bp 单元 × ~522 拷贝） | P42 报告只泛称 "repeat context" | 无法向客户解释"为什么这里的单 read junction 大多不可信" |
| B6 | **分析历史上曾设 60× sensitivity gate**（~15× 数据不可能达到，使所有候选"未过 gate"），交付版虽已去掉，但正文仍大量"not established / not an absence call"式否定表述 | 报告 §2 第 5 条、§6.6 | 结论只有"什么都不能说"，对客户无信息量 |

### 3.2 客户报告写法

| # | 缺点 | 证据 | 违反的规则 |
| :--- | :--- | :--- | :--- |
| B7 | **§2 没有"对你意味着什么"**，6 条全是方法学限定 | 报告 §2 | report skill：§2 Core Conclusions 是"逐对象该怎么办，大白话" |
| B8 | **内部审计术语密集**：QNAME、standard-eligible、comparator zero、interpretation boundary、raw-support/depth proxy | 报告 §4、§6.6、§8.2 | 客户是实验生物学家；工具字段名首次出现须给全称（QNAME 未定义） |
| B9 | **客户正文与图标题出现内部项目编号** P41/P42 | 报告标题"P42 whole-genome sequencing …"；Fig14 标题"P42 answers the P41 Section 7.3 customer questions" | 客户不知道 P41/P42 是什么；应称"the report dated 2026-10-09" |
| B10 | **header 缺 Prepared for** | 报告头 | report skill header 规范 |
| B11 | **样本表只列 5 个 DLD1 样本，但声称分析 18 个**；GM12878/Mutu I/DG75 无任何结果 | 报告 §3 | 读者无法知道其余 13 个样本得出了什么 |
| B12 | **同一限定语反复出现**（"not independent validation"、"repeat-context" 等在 §2、§4、§6、§8 多次） | 全文 | 用户规则"一篇文档内不做简单重复" |

### 3.3 图件与交付组织

| # | 缺点 | 证据 | 影响 |
| :--- | :--- | :--- | :--- |
| B13 | **Fig10 文字压在数据点上**；x 轴写 "chr11 anchor coordinate" 但右侧文字给的是 partner 坐标（chr21/chr3） | `figures/P42_Fig10_case_only_breakpoint_schematics.png` | 违反 figure skill §6 渲染后碰撞检查；坐标含义混淆 |
| B14 | **14 张图多为流程/审计示意**（方框箭头图 Fig13/Fig14），承载数据少 | Fig13、Fig14 | 图件数量多但信息密度低 |
| B15 | **846 GB 整基因组 BAM 放进报告目录**（硬链接，link count = 3，不额外占盘） | `genome_wide_data/bam/` | 报告包与原始数据传输混为一体；客户下载报告即被迫面对 846 GB；10-14 邮件已承诺 BAM 另走 Dropbox |

---

## 4. P44 自身的不足（对照后发现）

| # | 不足 | 实测 | 对应 P42 优点 |
| :--- | :--- | :--- | :---: |
| C1 | **交付文件 header 带本机路径**（必须修） | S04 locus BAM header 6 处；Manta VCF 4 处；chr11 SNV VCF 7 处（§1 命令） | A1 |
| C2 | 无机读的客户问题覆盖表 | 只有报告 §8.2 文字 | A2 |
| C3 | 关键事件逐 read 证据没进交付 | `analysis/20261024_read_rescue/stage6_genotype/read_support.tsv.gz` 只在内部 | A7 |
| C4 | §8.4 建议 PCR，但没给引物 | 只给 `sequences/junction_and_reference_templates.fa` | A6 |
| C5 | 无逐文件数据清单与外部输入缺口表；"SNV 只做 chr11"只在 §8.3 一句话 | — | A3、A4、A5 |
| C6 | §9 没写整基因组 BAM 的去向 | 客户 Deliverables 清单包含 BAM 相关文件；10-14 邮件承诺 Dropbox | A5（以"写明去向"而非复制 BAM 的方式） |
| C7 | 抽取窗口 ±5 kb，P42 的 C1（chr11:114,564,286）与 C4（chr11:114,598,286）落在窗口外 | 运行记录 §4 | 范围差异，非缺陷；客户问的是 target |

---

## 5. 明确不借鉴的做法

| P42 做法 | 不借鉴的理由 |
| :--- | :--- |
| 先由 caller 候选封闭输入，再审少数候选 | 违反 §2A；本项目的关键发现恰恰不在 caller 输出里（P44J0332 caller = neither） |
| 客户正文以否定式限定为主 | 限定应放在 §8.3 Limitations 一处写全，正文给结论（用户规则：一篇文档内不简单重复） |
| 客户文字中使用 P41/P42 内部编号 | 客户不认识；用日期称呼前一份报告 |
| 把整基因组 BAM 放进报告目录 | 报告包与数据传输分离；§9 写明另行传输即可 |
| "standard-eligible / comparator zero"等自造术语 | 客户读不懂；P44 用 present / absent / underpowered 三态并在 §5 一次性定义 |

---

## 6. 改进计划（P44 交付包，就地重建）

**前提**：交付包生成于 2026-10-25、尚未经用户审阅，属"审阅前修订"，在原目录 `custom_research_report_supplementary_20261025/` 上由 `12_build_delivery.sh 2026-10-25` 重建，不新开目录；旧版因从未被审阅、且可由同一脚本完全重现，不另归档（若用户要求对照，再从 git `9df28f8` 取出报告源重建）。

| 优先级 | 编号 | 改进项 | 具体做法 | 新增/修改文件 | 验收判据 | 工作量 |
| :--- | :---: | :---: | :--- | :---: | :---: | :---: |
| **必须** | P1 | header 去本机路径 + 审计 | 新脚本 `scripts/12a_sanitize_headers.py`：locus BAM 用 `samtools reheader` 重写 `@PG CL`（保留 ID/PN/VN，路径替换为文件名）；VCF/VCF.gz 用 `bcftools reheader` 去掉含路径的 `##` 行；输出 `methods/header_sanitization_audit.tsv`（文件、去掉行数、record 数前/后、body SHA-256 前/后） | `12a_sanitize_headers.py`；`12_build_delivery.sh` 调用它；`13_self_audit.py` 新断言 | 交付目录内 `grep -c '/home/\|/mnt/\|/Work_bio'` 对全部 BAM/VCF header 为 0；record 数与 body hash 前后完全一致 | 小 |
| 高 | P2 | 客户问题覆盖表 | `tables/client_request_coverage.tsv`：客户 4 问 + Analysis 2 的 3 个支撑项 + Deliverables 清单逐项 → 答案状态 → 证据文件 → 边界；报告 §1 末尾加同一张表（§8.2 只指向它，不复制） | `10_build_tables.py`；`_client_report_source.md` | 客户 10-13 邮件中每一个问句与每一项 Deliverable 在表中各占一行；self-audit 断言表中每个证据文件存在 | 小 |
| 高 | P3 | 关键事件逐 read 证据账本 | 对 Section 7.4 的 4 个 event 输出每条支持 read end：样本、read 名、read end、原始 BAM 中各 record 的 FLAG/坐标/CIGAR/MAPQ/SA/XA、模板 AS、基因组 AS、是否 duplicate | `tables/key_junction_read_evidence.tsv` | 行数 = 4 个 event 的支持 read end 合计（7 + 4 + 6 + 1 = 18，以 Stage 6 实际为准）；与 Fig5 绘制的 read 数逐 event 相等 | 小 |
| 高 | P4 | in-silico 引物 | Primer3 在 junction 模板上设计跨 junction 引物对；`bwa-mem2 mem -a` 与 BLAST 对 GRCh38 做全长精确命中计数，要求每条引物唯一命中预期位置；引物须避开 42 bp repeat 单元（否则多位点命中）；P44J0332 另给"chr11 侧引物 + (TTAGGG)n 引物"的 STELA 类方案说明 | 新脚本 `scripts/14_primer_design.py`；`tables/primer_design.tsv` | 每对引物两端唯一命中；产物长度在模板上可计算；repeat 内无法唯一的明确标注"无唯一引物" | 中（Primer3 若未装，按常设授权先列名称/用途再装） |
| 中 | P5 | 数据清单与外部输入缺口 | `tables/data_inventory.tsv`（交付内相对路径、bytes、md5、索引）；`tables/external_input_status.tsv`（junction PCR、FISH、全基因组 SNV/indel、整基因组 BAM：状态 + 对解读的影响） | 由 `12_build_delivery.sh` 生成 | 清单行数 = 交付目录文件数；md5 与 `md5sums_delivery.txt` 一致 | 小 |
| 中 | P6 | §9 写明整基因组 BAM 去向 | §9 表末加一行："Whole-genome BAM files of all 18 samples: provided separately via the shared transfer folder described in our 2026-10-14 message (about 850 GB)"——**措辞需用户确认传输状态后再定** | `_client_report_source.md` | 用户确认 | 极小 |
| 可选 | P7 | 上下文窗口扩到 ±50 kb | 仅当用户希望与 P41 窗口完全对齐时：W 外再加 context window，重跑 Stage 2–8（约 2.5 小时机时），C1/C4 类事件进 `junction_table_window_context.tsv` | `0_common.sh` 新参数；新日期分析目录 | 需用户决定（改变 data universe，§2A 要求开跑前确认） | 中 |

**执行顺序**：P1 → P3 → P2 → P5 → P4 → P6（待确认）→ 重建交付 → `13_self_audit.py`（新增 P1/P2/P3/P5 断言）→ `check_report.py` → 人眼看新增表头与报告改动处。

**报告改动范围**：只增加 §1 覆盖表、§8.4 引用引物表、§9 新文件行与 BAM 去向；科学结论、数字、分级不变（self-audit 保证）。

---

## 7. 需要用户决定的事

| # | 问题 | 默认 |
| :--- | :--- | :---: |
| D1 | 是否执行 P1–P5 | 执行（P1 为必须） |
| D2 | 整基因组 BAM 当前传输状态（Dropbox 是否已同步完成），用于 P6 措辞 | 待用户告知 |
| D3 | 是否做 P7（±50 kb） | 不做 |

---

## 8. 对今后项目的通用教训（拟沉淀到 delivery skill，待用户同意）

1. **交付前 header 去路径应成为固定步骤**：凡交付 BAM/VCF，交付脚本必须去 `@PG CL` / `##commandline` 等含本机路径的行，并以 record 数 + body hash 证明数据未变。可把检查加入 `delivery` skill 的"交付必给清单"与 self-audit。
2. **"客户问题覆盖表"作为每个定制分析交付的标配机读文件**，与报告 §8.2 互为索引。
3. **建议做 PCR 时随附 in-silico 引物与特异性筛查**，否则建议不可执行。
4. **审计严谨 ≠ 结论空洞**：限定语集中在 Limitations 一处；正文先给可行动的结论。P42 的教训是两者被混为一谈。
5. **候选封闭式流程即使后续审计再严，也补不回被排除的证据**：本项目 P42 与 P44 的结论差异，根源就在 data universe 的定义。

---

## 更新记录

- 2026-10-25 — 创建。
