# P45 chr11 locus junction 分析：P42 + P44 合并版（Vedita 18 WGS）— 研究计划与审阅

- **Plan Date（创建）：** 2026-10-25（文件名日期 = 创建日，永不改名；修订记在下方更新记录）
- **状态：** **计划待批准。** 本轮没有启动关键词“自动 / 全流程 / 全过程”，所以只建了目录、样本表和本计划。没有复制代码或数据，也没有运行任何分析。
- **样本表（SSOT）：** [`scripts/0_sample_table.csv`](../scripts/0_sample_table.csv)（从 P44 样本表复制，列名 `P44_role` 改为 `P45_role`，其余逐字相同）
- **客户要求原文：** Xin 2026-10-13 邮件，存于 `/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc/docs/P41_to_42_request.md`（执行阶段复制到本项目 `docs/`）
- **本计划依据的两份对比复盘（不在本文重复，只按节号引用）：**
  - 【P44←P42】`/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc/docs/P44-learned-from-P42/P44_learned_from_P42_1025.md`
  - 【P42←P44】`/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis/docs/P42-learned-from-P44/P42_learned_from_P44_20261025_CN.md`

## 更新记录

- 2026-10-25 — 创建。
- 2026-10-25 12:05 — **用户批准全部推荐**（Q1 = a 不设 60× gate；Q2 采用 structurally coherent 级；Q3 = b，S01–S10 genome-wide SNV/indel 联合 calling；Q4 = b，加做 G2；Q5 = a，单一完整包），并授权**全流程**。另授权：若额度用完，由定时任务在窗口后自动恢复本流程（登记在 `Agent_Assets/runtime/claude_code_active_project_autostart.json`，重试规则为重置时间 +5 min，兜底 3/4/5 h）。G1 不依赖 P44，先于 Stage 0 启动。
- 2026-10-25 16:00 — **§4.2 补一级（判据空档，不是改判据）**：§4.2 的 claim ladder 没有为 P44 的 Well_supported 类定义层级。第 2 轮数据中出现 2 个 Well_supported 事件（J1365、J0533：W 侧翼 chr11:114.5908 Mb，GM12878 全部 5 个样本均有，8–10 fragments）。9x v1 把它们并入 limited，低估了证据强度。现新增 `well_supported_junction_candidate`（措辞为 well-supported junction candidate, not independently confirmed），其余规则不变。v1 产物原样归档于 `superseded_for_comparison/analysis/v1_20261025_claim_ladder_well_supported_folded_into_limited/`。

- 2026-10-25 16:35 — **全流程完成**：交付包 `custom_research_report_20261025/`（212 个文件，2.0 G）已生成，`13_self_audit.py` 与 `check_report.py` 均通过；**待用户审阅，未发送**。运行中的发现与修正（SV 记录 END 丢失，影响 P44 Stage 8；G2a2 脚本 bug；Fig5 虚线与图例重叠）见运行记录 §3。

---

## 0. 结论先行

1. **代码以 P44 为骨架。** P44 的 data universe 定义正确：从 locus window 出发，经 mate、SA、XA 和双端 unmapped 五条路径 rescue，再让 junction template 与全基因组最佳 placement 竞争。P44 找到了 P42 漏掉的两个关键事件：S04 的 chr11→(TTAGGG)n 和 S05 的 376 bp fold-back。两份复盘都认可这一点（【P44←P42】§3.1 B1–B2；【P42←P44】§3.1–3.3）。
2. **P42 的贡献在三处：** 交付工程（header 去路径审计、问题覆盖表、数据清单、外部输入缺口表、逐 read 证据账本、引物特异性筛查）；多出的两个 caller，即 GRIDSS（18 样本联合）和 SvABA（DLD1）；结论分层的克制，也就是 claim ladder，同一数据内的一致性不写成独立验证。
3. **不需要改的科学计算一律复用，不重跑。** P44 第 2 轮（W = target ±50 kb，`analysis/20261025_read_rescue/`）此刻正在 tmux `p44r2_scan` / `p44r2_chain` 里跑。截至 2026-10-25 11:54，Stage 2 已扫到第 17 个样本，Stage 3 也在进行。P45 等它 `CHAIN_DONE` 后复制 Stage 1–8 的产物，再复制 P44 第 1 轮与 W 无关的 Stage 9 产物。P42 的 GRIDSS/SvABA 标准化结果和 C1–C4 结果只读复制，用作注释和对照。
4. **P45 新做的部分**不改变 P44 的 Stage 1–7 判据，只在其上增加五层：4-caller concordance、claim ladder 与 structural-coherence 层、逐 read 证据账本、引物设计与特异性筛查、合并版交付。
5. **开跑前需要你拍板的 5 个问题见 §8。** 其中 Q1、Q2 改变结论层级或 data universe，按 §2A 必须事先确认。

---

## 1. 背景与范围

### 1.1 三个项目的关系

| 项目 | 执行者 | 核心做法 | 交付状态（2026-10-25） |
| :--- | :---: | :--- | :---: |
| P41 `41_Vedict_18_wgs` | — | 18 样本全基因组 SV 初报（报告日期 2026-10-09，客户已收到） | 已发出 |
| P42 `42_Vedita_18_wgs_advanced_analysis` | Codex | full-window inventory → 282 candidate observations → 固定 C1–C4 → 严格 gate | `custom_research_report_20261024/` 已冻结，**未发送** |
| P44 `44_wgs_advanced_analysis_from_P41_by_cc` | Claude Code | 5 路 rescue → 全 locus discovery → 竞争性 re-genotyping → 分级 + depth/BAF | 第 1 轮包未发送；第 2 轮（±50 kb + P1–P7）运行中 |
| **P45（本项目）** | Claude Code | P44 骨架 + P42 的交付工程、额外 caller 与 claim ladder | 计划待批准 |

P45 **不写入** P41/P42/P44 的任何目录，只读它们的产物。

### 1.2 客户问题（Xin 2026-10-13），作为完成标准

| # | 客户原话要点 | P45 的回答载体 |
| :--- | :--- | :---: |
| A1-a | base-pair resolution breakpoints, microhomology / inserted sequence | `junction_table_target.tsv` + 逐 read 账本 |
| A1-b | junction-spanning read counts per event, before and after rescue | `junction_by_sample_long.tsv`（fragment 与 read-end 两列并存） |
| A1-c | well supported vs ambiguous because of repeat content | P44 分级 + claim ladder（§4） |
| A1-d | Section 7.3 events called by both callers vs one | 4-caller concordance（§5 Stage 8x） |
| A2 | +Dox clones vs S01–S03, only in +Dox; criteria（tier, read count, AF threshold） | per-sample state + detection power（§4.3） |
| D1 | SV VCFs from Manta and DELLY, unfiltered and filtered | 已有 VCF，去路径后交付（§5 Stage 12） |
| D2 | SNV/indel VCFs | **见 §8 Q3**：目前两个项目都只有 chr11 |
| D3 | annotated per-event TSV with caller column, regenerated for the new comparisons | **见 §8 Q4**：locus 级已有；genome-wide 版两个项目都没做 |
| D4 | transfer method and approximate size | 整基因组 BAM 已交付且客户确认收到（P44 D2）；§9 写明即可 |

---

## 2. 输入与数据量

### 2.1 原始输入（全部只读）

| 输入 | 路径 | 量（实测） |
| :--- | :---: | :---: |
| 18 个 BAM（bwa-mem2 2.2.1 `-Y`，MarkDuplicates 只标记不删除） | `/home/gao/projects_2026M7/41_Vedict_18_wgs/analysis/bam/S{01..18}_*.markdup.bam` | **18 个合计 845.12 GiB**（`stat -c %s` 逐个求和，2026-10-25 11:55）；S01–S10 481.99 GiB；DLD1 S01–S05 238.87 GiB；逐样本见样本表 `bam_GiB` 列 |
| Reference | `/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta`（与 BAM 同一份） | — |
| Manta / DELLY | `41_Vedict_18_wgs/analysis/sv/{manta_gw,manta_region_joint_20261007,delly}/` | 18 + 1 / 18 + 1 |

样本角色沿用 P44 的用户决定（2026-10-24 Q2）：S01–S10 逐样本报告；S11–S18 只作 artefact background panel，不逐样本报告。

### 2.2 复用的中间产物（执行阶段复制，均为只读来源）

| 来源 | 内容 | 体积 | 前提 |
| :--- | :--- | :---: | :---: |
| P44 `analysis/20261025_read_rescue/stage1–stage8` | W ±50 kb 的 inventory、genome scan、rescue BAM、clusters、assembly、genotyping、annotation、Manta/DELLY concordance | 约 4–5 GB（第 1 轮 ±5 kb 的 Stage 1–8 正式产物约 4.1 GB；第 2 轮截至 11:54 已 2.3 GB） | P44 第 2 轮 `CHAIN_DONE` 且 EXIT=0 |
| P44 `analysis/20261024_read_rescue/stage9_supporting` | 深度断崖 + bootstrap CI、chr11 SNV（S01–S10 联合）、BAF/LOH、VCF 整理 | 184 MB | 与 W 无关，P44 第 2 轮同样沿用 |
| P42 `analysis/20261023_fullcohort_repeat_resolve/07_external_sv_dld1_normalization_v1/` | GRIDSS 18 样本联合 + SvABA DLD1 germline/somatic 标准化事件表 + audit.json | 0.5 MB | 已完成 |
| P42 `custom_research_report_20261024/results/` 的 C1–C4 相关表与 `primer_specificity.tsv` | comparison annotation；C2 引物作对照 | < 10 MB | 已冻结 |
| P42 / P44 脚本 | 见 §6 复用清单 | — | 复制为快照，记录 md5 |

磁盘：`/mnt/ex_8T_SSD` 剩余 2.6 T，充足。复制时不使用硬链接，避免与源项目共享 inode 后被误改。

---

## 3. 证据模型 decision table（§2A；正式运行前须经审阅）

**Data universe 沿用 P44 第 2 轮，不再改动：**
- 抽取窗口 W = chr11:114,554,000–114,675,500（target ±50 kb，与 P41/P42 窗口一致，P42 的 C1/C4 落在其中）。
- Rescue 集合 = W 内 record ∪ mate 在 W ∪ SA 段在 W ∪ XA hit 在 W ∪ 双端 unmapped pair 经 W 序列筛选命中；对每个 QNAME 取回全部 record。
- 目标区 junction 的定义是：一端落在 target ±50 bp 内。

| 研究问题 | 证据单位 | 保留 / 去重 | 判定指标 | 允许的结论 | 不允许的结论 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| 位点有哪些 junction | record → read end → fragment | 沿用 P44（全部 record 进 audit；dup 保留但不计数） | 沿用 P44 Stage 4 | 候选 junction + 全部支持 record | 用 P41/P42 是否报告过来判断真假 |
| rescue 前后支持数 | fragment（主）+ read end（并列） | 沿用 P44 Stage 6 | before = P41 口径；after = 竞争性模板支持 | 变化量及来源路径 | 把 read end 数当分子数 |
| 断点 bp 结构 | split read 序列 + contig | 沿用 P44 Stage 7；**新增**逐 read query overlap/gap（P42 C2 的做法推广到全部 ≥2-fragment 事件） | 坐标、方向、microhomology、insert、不确定区间 | 断点 ± 区间 | 在区间内给单一坐标 |
| repeat 歧义 | flank 的全部 placement | 沿用 P44 | AS margin、拷贝数 | Ambiguous (repeat) | 以 primary 或 MAPQ 代替定位判断 |
| caller concordance | junction × sample × caller × record | **新增** GRIDSS（18 样本）、SvABA（DLD1）；Manta/DELLY 沿用 | 两端 ±100 bp 且方向一致 = concordant；单端 = anchor-only；另报最近距离，也报是否满足 P42 的 ±50 bp | 每个 caller 的状态 | caller 一致 = 独立验证 |
| clone-specific | fragment | 沿用 P44 Stage 6 统一基因分型 | absent / underpowered / present + 检测功效（§4.3） | +Dox-only / shared / underpowered | 低深度下的 0 = 不存在 |
| **结构一致性（新增，见 §4.2）** | junction 位置 vs 深度断崖 CI vs LOH 区段 | Stage 9 产物，预先写定判据 | 落入 / 不落入 | “与局部结构图式一致” | 独立生物学验证、Dox 因果 |

**Filter-impact audit：** P44 的 `evidence_filter_impact.tsv` 原样继承；新增层（caller 匹配、claim ladder 升降级）各出一张 before/after 表。

---

## 4. 判定层级（预注册；结果出来后不改，改则新版本并归档）

### 4.1 证据强度：沿用 P44 分级，不改动

Well supported / Limited / Ambiguous (repeat) / Weak / Not supported after rescue。定义见 P44 计划 §4 及其 2026-10-25 更新记录第 6 条，P45 不复制定义，报告 §5 一次性给出。

### 4.2 结论层级：新增 claim ladder（来自【P42←P44】§6.3）

这是 P44 分级之上的**另一列** `claim_level`，不替换 P44 分级。

| claim_level | 最低证据 | 报告措辞 | 禁用推论 |
| :--- | :---: | :---: | :--- |
| Mapping observation | 1 fragment，或只有 alignment/contig 线索 | weak observation | true locus、clone-specific |
| Repeat-ambiguous junction | 有支持，但一端或两端有多处等价 placement | junction-like evidence with ambiguous placement | 具体 partner 染色体 |
| Limited junction candidate | ≥2 fragments，完成模板/基因组竞争，报告了不确定区间 | limited junction candidate | 已确认、因果、完整重排 |
| Structurally coherent limited candidate | Limited 或 Ambiguous (partner end)，**并且**满足下面两条：① chr11 侧断点区间与本样本深度断崖 95% CI 有交集；② 断崖远端在同一样本为 LOH（≥90% baseline-het 位点变为 homozygous） | limited short-read evidence coherent with the local copy-number pattern | 独立验证、Dox 因果、完整构型 |
| Confirmed | 上一级 + 独立 assay（PCR/FISH 原始结果） | confirmed | — |

本数据没有任何独立 assay 输入，所以 **Confirmed 一级预期为 0**。这一点写进 `external_input_status.tsv`，不当作阴性结论。

### 4.3 P42 的 “60× sensitivity gate” 如何处理（§8 Q1）

P42 要求 ≥3 standard-eligible QNAME、局部深度 ≥60×、有 sequence-level junction，三项都满足才算 high-confidence。约 15× 的数据实际上不可能满足 60×，所以这个 gate 下所有候选都“未通过”，对客户没有信息量（【P44←P42】B6）。

**推荐：** 不作为分类 gate。由 P44 的逐样本检测功效替代：在 AF 0.1 / 0.25 / 0.5 下观察到 ≥2 和 ≥3 fragment 的 Poisson 概率。这让“未检出”可以读成“在该 AF 以上未检出”。P42 的 standard-eligible 计数作为一列保留（C1–C4 对照用）。

---

## 5. 分步流程

所有 P45 脚本放在 **`scripts/20261025/`**（日期 = 本轮运行日；按【P42←P44】§4.5，产出交付的脚本按日期冻结）。共享配置 `scripts/20261025/0_common.sh` 由 P44 的 `0_common.sh` 改写（PROJ 改为 P45，`P44_*` 环境变量改为 `P45_*`）。

| Stage | 脚本 | 内容 | 输出（`analysis/20261025_*`） | 机时 |
| :--- | :---: | :--- | :---: | :---: |
| 0 | `0_wait_and_import.sh` | 等 P44 第 2 轮 `CHAIN_DONE`；记录 P44 第 2 轮所跑脚本的 md5，并与 projects_2026Q4 git 中 P44 的最新提交比对（有未提交改动则在运行记录中写明差异）；复制 Stage 1–8（±50 kb）、Stage 9（20261024）、P42 07_ 与 C1–C4 表；写 `import_manifest.tsv`（源路径、源 md5、目标 md5、源脚本 md5） | `20261025_inherited/` | 等待约 2 h + 复制 < 10 min |
| 8x | `8x_caller_concordance_4callers.py` | 在 P44 `8_caller_concordance.py` 基础上增加 GRIDSS/SvABA；逐 junction × 样本 × caller 记录 record ID、匹配类别、最近距离、±50/±100 双标志；对 C1–C4 与 P42 的 0/20 结果互相校验 | `stage8x_concordance/` | < 30 min |
| 9x | `9x_structural_coherence.py` | 按 §4.2 判据把 junction 区间与深度 CI、LOH 区段做交叉；输出每事件 `claim_level` 与判定依据 | `stage9x_coherence/` | < 10 min |
| 10 | `10_build_tables.py`（P44 改） | 主表增加 4-caller 列、`claim_level`、P42 crosswalk 列（C1–C4、standard-eligible）；**新增** `client_request_coverage.tsv`、`external_input_status.tsv` | `stage10_tables/` | < 10 min |
| 10b | `10b_read_evidence_ledger.py` | 每个 ≥2-fragment 事件（预计 4 个 + C2 对应事件）的每条支持 read end：匿名 read ID、R1/R2、strand、原始与再比对 record 的 FLAG/坐标/CIGAR/MAPQ/SA/XA/SEQ、模板 AS 与基因组 AS、dup 状态、query overlap/gap 与 overlap 序列、fragment 分组 | `key_junction_read_evidence.tsv` | < 10 min |
| 10c | `10c_p42_reconciliation.py` | C2 / EV0623 / P44J0289 三方对照（坐标、计数单位、universe、判据），以及 C1、C3、C4 在 ±50 kb 下的 P45 状态；**内部表，不进客户包** | `docs/P45_P42_P44_reconciliation.tsv` | < 5 min |
| 11 | `11_figures.py`（P44 改） | 按【P42←P44】§4.3/4.4/4.6 修：Fig5 B 端标注改为 “telomeric-repeat sequence; representative placement shown” 并加 ambiguous partner 标签；junction view 加 read ID、R1/R2、strand、fragment 分组；新增窄版 client-summary panel；生成 `figure_manifest.json`；渲染后检查文字重叠（figure skill §6） | `stage11_figures/` | < 15 min |
| 14 | `14_primer_design.py`（P44）+ `14b_primer_specificity.py`（P42 `2_check_primer_specificity.py` 改） | 为每个 ≥2-fragment 事件和 C2 设计跨 junction 引物；`bwa-mem2 -a` + 全长精确命中计数要求唯一命中；避开 42 bp repeat 单元；telomeric junction 另附 STELA 类方案说明；无法唯一时明确写“无唯一引物” | `stage14_primers/` | < 30 min |
| 12a | `12a_sanitize_headers.py`（P44）+ DELLY BCF→VCF 审计（P42 做法） | BAM `@PG CL`、VCF `##commandline` 等去本机路径；审计 record 数前后一致、VCF body SHA-256 前后一致 | `methods/header_sanitization_audit.tsv`、`delly_conversion_audit.tsv` | < 20 min |
| 12 | `12_build_delivery.sh <YYYY-MM-DD>`（P44 改） | 报告源 `docs/_client_report_source.md` → `custom_research_report_YYYYMMDD/`；`data_inventory.tsv`（逐文件 bytes、md5、索引）；`md5sums_delivery.txt` | 交付目录 | < 15 min |
| 13 | `13_self_audit.py`（P44 改） | 原有断言 + 新断言：header 路径为 0；问题覆盖表中每个证据文件存在；账本行数 = 各事件 read end 合计 = 图中 read 数；caller 列与 8x 表一致；报告中每个数字从当前 TSV 重算；clientside 文字不含 P41/P42/P44/P45 编号；`check_report.py` 通过 | `self_audit.txt` | < 5 min |

**可选 Stage（按 §8 决定）：**
- **G1 genome-wide SNV/indel VCF（Q3）：** `bcftools mpileup -q 20 -Q 20 | call -mv`，按染色体并行，联合 calling。范围 S01–S10 或 18 样本。预计 6–10 h，线程 ≤ 56。产出 `.vcf.gz` + 去路径 header。
- **G2 genome-wide clone-specific per-event TSV（Q4）：** 基于 P41 已有的 Manta/DELLY genome-wide 调用，按客户 A2 口径（S04/S05 vs S01–S03）重新生成带 caller 列的事件表，并写明 tier、read 数、AF 阈值。预计 2–4 h，主要是写代码和审计。

**资源与安全：** 持续 ≤ 28 核 / 56 线程。P44 第 2 轮结束前，P45 不跑任何重计算。长任务在 tmux `p45_<stage>` 里跑。只写 P45 下的新日期目录。中途打断的后果：Stage 0 可重入（逐文件 md5 校验后跳过）；其余 stage 都在分钟级，重跑即可。

---

## 6. 代码复用清单

| 来源 | 文件 | P45 处理 |
| :--- | :--- | :---: |
| P44 | `p44lib.py`、`1_`–`8_`、`9b_`/`9c_`/`9d_`、`2_*.sh`、`3b_*.sh`、`0_run_rescue_chain.sh` | **原样复制为快照**（不运行，只为可追溯；这些 stage 的产物从 P44 复制） |
| P44 | `10_build_tables.py`、`11_figures.py`、`12_build_delivery.sh`、`12a_sanitize_headers.py`、`12b_delivery_tables.py`、`13_self_audit.py`、`14_primer_design.py`、`docs/_client_report_source.md` | 复制后修改（§5） |
| P42 | `20261024/2_check_primer_specificity.py`、`20261024/10_build_unified_client_package.py` 中的 header 脱敏与 DELLY 转换片段、`20261023/7_normalize_external_sv_and_dld1.py`（只读参考） | 摘取逻辑并改写进 P45 脚本；来源行号记入运行记录 |
| P42 | `results/customer_request_coverage.tsv`、`data_inventory.tsv`、`external_input_status.tsv` 的列结构 | 作模板，内容由 P45 重新生成 |

复制规则：快照放 `scripts/20261025/from_P44/`、`scripts/20261025/from_P42/`，文件名不变。被修改的版本放 `scripts/20261025/` 本层，可以直接 `diff` 对照来源。

---

## 7. 交付设计

- **单一交付包**：`custom_research_report_YYYYMMDD/`，完整回答 10-13 请求（§8 Q5 确认后定名）。报告全英文，按 §1–§9 结构，三处日期一致，署名为用户。
- **客户文字的写法：**
  - 不出现 P41/P42/P44/P45。前一份报告称作 “the report dated 2026-10-09”。
  - 限定语只在 §8.3 Limitations 写全一次，正文给结论（【P44←P42】B7、B8、B12）。
  - §2 按样本组写“结果 / 对你意味着什么”（P44 写法）。
  - §3 写 S11–S18 背景面板得出了什么，避免 P42 B11 的问题。
- **P42 已冻结报告里的结论**（C1–C4 未过 gate）与 P45 不同。按 report skill 的归因框架，写作“进一步分析在更完整的 rescue universe 中细化了证据”。该报告从未发出，客户不知道它存在，所以客户文字里不提及它。内部对照表 `docs/P45_P42_P44_reconciliation.tsv` 存档备查。
- **包内容：** 报告 + `tables/` + `figures/`（含 manifest）+ `locus_BAM/`（每样本 rescue BAM，已去路径）+ `sequences/`（模板、contig、引物）+ `supporting_analyses/` + `genome_wide_SV_VCF/`（Manta/DELLY filtered 与 unfiltered，DELLY 转成 VCF）+ `methods/`（审计表）+ md5。**整基因组 BAM 不放入包中**，§9 写明已另行交付。
- **只准备、不发送**；生成后交给用户过目，附逐条主张对数据来源的对照表。

---

## 8. 开跑前需要你拍板的问题

| # | 问题 | 选项 | 推荐默认 | 为什么需要你定 |
| :--- | :--- | :--- | :---: | :--- |
| Q1 | P42 的 60× high-confidence gate | (a) 不作 gate，用检测功效替代，P42 计数作为列保留；(b) 作为 claim ladder 之外的独立一列，同时报告 | **(a)** | 改变结论层级（§2A） |
| Q2 | claim ladder 中 “structurally coherent” 一级的判据（§4.2 ①②）是否采用 | 采用 / 只报告事实不设该级 | **采用** | 新增结论层级，必须事先预注册 |
| Q3 | 客户 D2 “SNV/indel VCFs”（genome-wide） | (a) 不做，§9 和缺口表写明只提供 chr11；(b) S01–S10 genome-wide 联合 calling（约 6–10 h）；(c) 18 样本（约 10–16 h） | **(b)** | 交付范围；改变总耗时 |
| Q4 | 客户 A2 / D3 的 genome-wide clone-specific per-event TSV | (a) 只做 locus 级（与 P44 相同）；(b) 加做 G2 genome-wide 版 | **(b)** | 客户原文要的是 genome-wide；两个项目都没做 |
| Q5 | 交付形态 | (a) 单一完整包 `custom_research_report_YYYYMMDD/`，取代未发出的 P42/P44 包；(b) supplementary 包 | **(a)** | 决定报告结构与命名 |

另有两件事按默认做，不需要拍板（改一行即可回改）：
- caller 匹配以 ±100 bp + 方向一致为主口径，另附 ±50 bp 标志，以便与 P42 对得上。
- 若 P44 第 2 轮失败或在今晚之前未完成，P45 不等待：以快照代码在 P45 内自行跑 Stage 1–8（约 3 h），并在运行记录中注明。

---

## 9. 完成定义

1. 交付目录已生成；`13_self_audit.py` 全部断言通过；`check_report.py` 通过。
2. §1.2 表中每一行在 `client_request_coverage.tsv` 中都有对应行，状态为“已回答”或“不可得 + 原因”。
3. 运行记录 `docs/P45_key_operations_and_parameters_log_1025.md` 写全：每条决定配 file:line 和独立复核命令。
4. `Agent_Assets/runtime/active_projects.json` 中状态改为 `complete`。
5. 交给你过目，附主张对证据的对照表。**不对外发送。**

## 10. 预计时间（批准且 P44 第 2 轮按时结束的情况下）

| 段 | 墙钟 |
| :--- | :---: |
| 等待 P44 第 2 轮 + 导入 | 约 2 h（预计 13:30–14:30 完成） |
| 写 / 改脚本 + Stage 8x–14 | 3–4 h |
| G1（若选 b） | 6–10 h，可与上一段并行 |
| G2（若选 b） | 2–4 h |
| 交付构建 + 自审 + 报告 | 2–3 h |
