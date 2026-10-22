# P42 核心运行记录与参数选择

- **创建日期**：2026-10-14（文件名日期即创建日，永不改名；后续修订在文末「更新记录」追加）
- **项目**：`/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis`
- **前置项目**：项目 41，`/home/gao/projects_2026M7/41_Vedict_18_wgs`（= `/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs`）；已交付 `custom_research_report_20261009/Chr11_Repeat_Region_WGS_SV_Report_1009.md`，**冻结，不动**
- **本文性质**：内部文档，只记录**已发生**的操作与决定；每条配出处与可独立复核的命令。估算与实测分开标注。
- **状态（2026-10-14）**：客户已要求继续深入研究，但**选项 A / B / C 尚未选定、预算尚未追加批准，分析尚未启动**。

---

## 1. 客户请求与回信（2026-10-13）

| # | 决定 / 事实 | 类别 | 出处与复核 |
| :--- | :--- | :---: | :--- |
| 1.1 | 客户（Xin）要求两项：**Analysis 1** 为 chr11:114,604,000–114,625,500 的 junction 级分析 + read rescue；**Analysis 2** 为 S04、S05（+Dox 克隆）对 S01、S02、S03 的克隆特异事件。并要求全基因组 BAM/CRAM、SV VCF（Manta + DELLY，过滤/未过滤）、SNV/indel VCF、带 caller 列的 per-event TSV；预算上限 **$840 = 7 h**（$120/h） | 事实 | 客户来信（用户贴入对话）；费率见 `/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs/docs/custom_request.txt` 中 "7 additional hours @ $120/hour = $840" |
| 1.2 | **估时 25–45 h（$3,000–5,400）**。拆分：Analysis 1 约 12–22 h，Analysis 2 约 6–10 h，交付物与报告约 7–13 h | **估算，非实测** | 用户定区间；拆分是我推出的。邮件脚本 `/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs/scripts/22_draft_reply_xin_1013_v12_transfer_user_wording.py` 的 BODY |
| 1.3 | **选项 B = Analysis 2 only，含 SV 与 SNV/indel VCF，12–20 h（$1,440–2,400）**；选项 C = Analysis 1 限于客户挑的 7.3 事件，固定约 6–8 h + 每个 junction 约 0.5–1 h（5 个约 9–13 h） | **估算，非实测** | 同上；我自己的独立估算 B 为 14–24 h（含 SNV/indel 人工 2–4 h，机器数天），邮件取用户定的 12–20 h |
| 1.4 | 此前 Option 2 报价里的 "3–4 h" 指标准流程调出事件的 **IGV 目视复核**，不含 read rescue、局部重组装；邮件中已写明这一口径差别，未否定旧估计 | 事实 | `docs/custom_request.txt`（项目 41）"the extra 3-4 hours for the highly noisy repeat regions"；邮件脚本 22 |
| 1.5 | 邮件**只写 Part 1（报告包，预计数百 MB，未实测）/ Part 2（全基因组数据）**，不写具体目录名，因目录日期要等用户审阅定稿 | 决定 | 脚本 22 的 "Deliverables and size" 段 |
| 1.6 | 数据传输写法：标准方式共享 Dropbox 文件夹，备选 SFTP；若单位屏蔽我方 SFTP 的服务器名与端口（已有哈佛系医院客户遇到过），可在家下载到外置硬盘。"terminal" 我按"端口"理解 | 用户措辞 + 我的理解 | 脚本 22；用户原句见对话 2026-10-13 |

## 2. 项目 41 的数据基线（本项目起点）

以下全部在 `/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs` 下实测（2026-10-13）。

| # | 事实 | 复核命令 |
| :--- | :--- | :--- |
| 2.1 | 全基因组 markdup BAM：18 个，合计 **846 GB**（`*.sorted.bam` 另有 706 GB，不属交付候选） | `du -ch analysis/bam/*.markdup.bam \| tail -1`；`ls analysis/bam/*.markdup.bam \| wc -l` |
| 2.2 | 已交付的只有 locus 窗口子集 BAM：36 个文件（18 对 BAM + bai），共 33 MB；**全基因组 BAM 尚未交给客户** | `ls custom_research_report_20261009/region_bam \| wc -l`；`du -sh custom_research_report_20261009/region_bam` |
| 2.3 | 区域 VCF 各只有 **1 条记录**，是同一个 chr11:114598292 的 **36 bp 插入**，FILTER=PASS，Manta 与 DELLY 的 18 个样本基因型完全一致：S01–S05（DLD1）0/1，S06–S10（GM12878）0/0，S11–S18（Mutu I、DG75）1/1 | `bcftools view -H custom_research_report_20261009/region_vcf/diploidSV.vcf.gz`；`bcftools query -f '%CHROM:%POS %SVTYPE %SVLEN [ %GT]\n' ...` |
| 2.4 | 该插入序列（正链，位于参考碱基 G 之后）`TGGGGTTCAGTCCCTGCAGCTGCTCTCAAGGGCTGA`；Manta 标注 `HOMLEN=2, HOMSEQ=TG`；距目标区域上游 5.7 kb；7.3 节对应 EV0572（chr11:114,598,293–294，与 VCF 记录差 1 bp） | 报告 §7.3、§7.5 与 `Chr11_Repeat_Region_WGS_SV_Report_1009.md` 第 69、236、283 行 |
| 2.5 | **未核实**：该插入是否为已知多态（dbSNP、gnomAD），是否落在重复序列或基因内 | — |
| 2.6 | 7.3 节窗口为 ±50 kb，共 4,266 个 junction 事件；tier A = 2，B = 37，C = 4,227；69% 标为 repeat-ambiguous；55% 的事件其伙伴端落在 chr2:32,916,000–32,917,000（视为 mapping artefact） | 报告第 225–231 行 |
| 2.7 | 目标区域内有 11 个 tier B 事件，各只在 1 个样本里，没有 tier A 支持 | 报告 §7.3 第二张表 |
| 2.8 | `Region_events_all`（4,266 行 × 32 列）**没有 caller 列**；"caller" 在报告里是 tier 判定用的内部标记（Manta 或 DELLY 在两端各 100 bp 内同样本调用）。**新的汇总表必须新增 caller 列** | `P41_SV_results.xlsx` 的 `Region_events_all` 列名；报告第 107 行 |
| 2.9 | **caller 一致性（更正我草稿里的一处错误）**：整个 ±50 kb 窗口内只有 EV0572 被 Manta 与 DELLY 同时支持；其余事件（含目标区域内 11 个 tier B）无 caller 调用。我最初把报告 2c 条目（仅覆盖 chr11:114.603–114.611 Mb）误当成整个 7.3 窗口，草稿已改正 | 脚本 `20_draft_reply_xin_1013_v10_caller_concordance_fixed.py` |
| 2.10 | **未逐行核对** `Region_events_all` 每个事件的 tier 判定里 caller 的取值，只核对了 VCF 层面 | — |
| 2.11 | 另有一个 `/home/gao/projects_2026Q4/41_Vedita_18_wgs`，只含 `analysis/`（含 `bam/` 和 `P41_SV_results_draft.xlsx`），我没有分析它，与上面的项目 41 的关系**未核实** | `ls /home/gao/projects_2026Q4/41_Vedita_18_wgs` |

## 3. 环境与工具

| # | 事实 | 复核命令 |
| :--- | :--- | :--- |
| 3.1 | 2026-10-13 经用户明确同意，在 `regular_bioinfo` 里装了 **openpyxl 3.1.5**（此前该环境没有；pandas 2.3.3） | `grep openpyxl /Work_bio/gao/configs/.conda/envs/regular_bioinfo/conda-meta/history` |
| 3.2 | ⚠ 本机 `conda run -n regular_bioinfo python` 实际落到 `/mnt/ex_8T_SSD/conda_envs/sc_tut_spatial`（Python 3.11.16），没有 openpyxl。**用环境里 Python 的绝对路径**：`/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python`（3.10.20）。bcftools 等同样用该环境 `bin/` | `conda run -n regular_bioinfo python -V` 与 `/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python -V` 对比 |
| 3.3 | 修改了邮件检查器 `/home/gao/.claude/hooks/checks/check_client_email.py`：脚本含 `# transfer-answer-intended` 时，只豁免 `dropbox`、`upload` 两个禁用词，其余判据不变。用户 2026-10-13 明确要求。该修改**尚未提交**；`test_checks.py` 当时 65/65 通过，但**没有为这条豁免新增用例** | `cd /home/gao/.claude && git diff hooks/checks/check_client_email.py` |

## 4. 邮件草稿的已知坑

- 用 IMAP APPEND 追加的草稿**没有 Date 头**；在 Gmail 里打开并保存后，Gmail 会重写正文：折行、把 HTML 加粗转成 `*…*`、加上 Date 头，并**改变草稿顺序**。2026-10-13 因此出现过用户看到旧写法的情况（旧版本被 Gmail 重新保存后留在草稿箱）。
- **最终版脚本**：`/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs/scripts/22_draft_reply_xin_1013_v12_transfer_user_wording.py`（脚本 11–22 为 12 个迭代版本；每次运行新增一封草稿，不删旧的）。
- 回读判据：用 IMAP `BODY.PEEK` 读 `[Gmail]/Drafts`，检查关键句是否在；不要靠脚本退出码判断。
- 本轮**邮件尚未发送**，是否发出由用户本人决定。

## 5. 待客户答复 / 待用户决定

1. 选项 A / B / C 选哪个；是否追加预算（三个选项都超过 $840）。
2. 是否传全基因组 BAM（约 846 GB）及怎么传（见 1.6）。
3. 启动正式分析需要用户给出"自动 / 全流程 / 全过程"之一；在此之前只做背景调研和研究计划。
4. 本项目还缺：`scripts/0_sample_table.csv`（样本表，作为 SSOT）、`P42_research_plan_and_review_MMDD.md`（研究计划）。

## 更新记录

- 2026-10-14 — 创建本文，补入项目 41 回信阶段（2026-10-13）的操作、核对与已知坑。
- 2026-10-22 — 用户以“自主完成P42分析全过程”授权正式启动 Analysis 1，并最终将范围定为 S01–S10（显式包含 S02、S03、S10；排除 Mutu I 与 DG75）。新建的 `docs/P42_research_plan_and_review_1022.md` 与 `scripts/0_sample_table.csv` 是本轮运行的计划与样本 SSOT；本文件保留 10-14 的历史记录，不将旧的“尚未启动”状态改写为当时已启动。
- 2026-10-22 — 正式 read-rescue 完成：P41 的 target-evidence 表按“至少 1 条 non-duplicate、两段 MAPQ >=20 的 split read”规则得到 179 个候选（S01–S10）；每个候选均重新提取到至少一个原始支持 fragment，并以 `bwa-mem2 mem -a -Y` 独立重比对。结果为 4 个 multi-fragment supported、173 个 single-fragment supported、2 个 not reconstructed、0 个 ambiguous/not evaluable。正式输出：`analysis/20261022/results/`；准备供用户审阅的补充交付：`custom_research_report_supplementary_20261022/`，尚未发送。
