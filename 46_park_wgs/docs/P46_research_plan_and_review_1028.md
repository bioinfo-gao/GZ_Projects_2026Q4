# P46 — Yeji Park 第二批 WGS（10 样本）：表单解读、待问清单与工时估算

- **Plan Date（创建，不可更改）**：2026-10-28
- **项目目录**：`/home/gao/projects_2026Q4/46_park_wgs/`
- **客户**：Yeji Park，HHMI / Boston Children's Hospital（TaekJip Ha lab）
- **前序项目**：P31 `/home/gao/projects_2026M7/31_Yeji_Park_RPE1_CRISPR_dual_cut_wgs/`
  （hTERT RPE-1，PPP1R2 双切 202 bp 单等位缺失 + 全基因组脱靶 + 载体整合筛查）
- **客户表单**：`Sample_Registration_Form_Yeji Park.xlsx`（项目根目录）
- **样本表（SSOT）**：`scripts/0_sample_table.csv`（本文不重复样本清单）
- **当前状态（2026-11-06 21:50 复核）**：客户已交回问卷（`client_materials/wgs_questionary.docx`，2026-11-06 13:13
  邮件附件），**§3 的 4 项待问已有 3 项得到答复**；**测序数据仍未到**（`/home/gao/Dropbox/` 下无
  `Park_Yeji`/`Park, Yeji` 目录，`find -L /home/gao/Dropbox -maxdepth 2 -iname '*park*'` 返回 0 命中）；
  深度已定 **45× raw**（客户 10-30 来信 "Please run 45X!" + Quote 10282602）；**未开跑**，本轮无授权关键词
  （"自动"/"全流程"/"全过程"），故本文只到计划层。

## 更新记录

- 2026-10-28 — 建档：表单解读、与 P31 的关系、待问客户 4 项、P31 工时实况、P46 工时估算。
- 2026-10-28 — 新增 §3.1 测序深度核查；§2 依据 1 改为引用 09-01 **实际发出版**（初稿误引脚本草稿）。
- 2026-10-29 — 用户出具 **Quote 10282602（WGS 45x × 10，$5,760）**，存 `client_materials/`；
  按 09-01 回信起草致客户邮件建议 45–50× raw，存 Gmail 草稿（未发送、附件待用户手动挂），
  存档 `docs/P46_email_depth_quote_1029.md`。另据 Kuo-Chan 10-27 来信：样本 10-20 前已到，
  本周建库、下周测序、预计 11-09 当周交付数据；Yeji 10-20 称 "It is the same samples before"（含义未确认）。
- 2026-10-29 — **用户已发出**（14:52 UTC，Sent `1a188259fa838bdc`），并加了一句**价格承诺**：
  靶点缺失、全基因组脱靶、CRISPR 造成的其他变化三项分析对 10 个样本**免费**。⇒ §5 中这三项对应的
  定制工时（约 15–22 h）不再向客户计费；载体筛查等未在信中承诺。原文与范围划界见
  `docs/P46_email_depth_quote_1029.md`「实发版」一节。

- 2026-11-06 — 客户交回问卷 + 两封邮件（PI = Taekjip Ha；"controls are same cell line"；"more single
  allele samples"；"both allele editing samples"）。新增 **§7 问卷解读与样本分组**、**§8 分析步骤设计**、
  **§9 待确认事项（已收窄为 4 项）**；`scripts/0_sample_table.csv` 补入 `client_global_no`、
  `expected_target_genotype`、`expected_del_alleles`、`prior_evidence`、`cell_line` 五列并改写 role 与 notes。
  §3 的四项待问：第 1、2、3 项**已答**（见 §7），第 4 项（浓度）仍未答。
- 2026-11-06 — **用户裁定关闭浓度问题**：Hung 已完成建库、测序即将出结果，该歧义不再影响任何可做的
  决定 ⇒ §3 第 4 项与 §9 第 3 项同时标记为已关闭，**不转 Jing、不问客户**。待确认事项由 4 项降为 3 项
  （Cas9 递送方式、样本编号对齐、是否与 P31 合并分析）。

---

## 1. 表单读到了什么（实测）

读取方式：P31 的富文本逐 run 解析脚本
`/home/gao/projects_2026M7/31_Yeji_Park_RPE1_CRISPR_dual_cut_wgs/scripts/1_inspect_registration_form.py`
（读深红色 FFC00000 作答）+ 全单元格转储（读普通文本）。

| 字段 | 内容 | 判定 |
| :--- | :---: | :--- |
| 样本 | 10 个 gDNA：Ctr_1–3、Sample 1–7 | 已填（详见样本表） |
| Service / Species | 全部 WGS（Whole-Genome Sequencing）/ Human | 已填 |
| Paired (B18) / Biosafety (B20) / Human-derived (B21) | 均以深红色标注 **No** | 已填，作答方式同 P31 |
| Sample Note (B16) | 模板原文 | **未填** |
| Analysis Type (B17) | 模板原文 | **未填** |
| Which Samples Are Control (B19) | 模板原文 | **未填** |
| Analysis Purpose (B22) | 模板原文 | **未填** |
| Thresholds (B23, 可选) | 模板原文 | 未填（按我方默认） |
| Quote Number | QTE-26-07-30-001 | **沿用 P31 旧单号** |
| Project / Expected Delivery | Dual Cut Project / Aug 10 2026 | **沿用 P31 旧表** |
| 浓度 (D42–D51) | `total 200ng/ul`，体积 20 µL | 写法有歧义，但**已不追问**（建库已完成）——见 §9 第 3 项 |

⇒ 表单回答了「送什么」，**没有回答「要做什么」**。决定范围与报价的 B16 / B17 / B19 / B22 全空。

## 2. 客户要做什么 —— 推断（未经客户确认）

- **推断**：P31「Dual Cut Project」的扩大版 —— 多个处理克隆 + 3 个对照，目标为
  PPP1R2 靶点（GRCh38 chr3:195,542,947 / 195,543,149）202 bp 单等位缺失，加全基因组脱靶。
- **依据**：
  1. 我方 2026-09-01 18:42 UTC 实际发出的回信（Gmail Sent；与 P31
     `scripts/38_shorten_client_reply_0901_depth.py` 草稿不同，用户发前改过）建议下一批
     "about 30x（~20x at least）usable … roughly 45x to 50x（~30x）of raw output"，并加未处理亲本对照；
     本批恰好出现 3 个 Ctr + 7 个处理样本。
  2. P31 `docs/with_client_info_20260825/P31_round2_plan_and_review_0825.md` §6.5 已记录
     「下一批（10 样本 / 4 笔业务）」—— 样本数一致。
- **未查到的**：Gmail 中 Yeji 最后来信 2026-08-31（仅致谢），我方最后去信 2026-09-01；
  此后无关于本批的往来（2026-10-28 回读）。⇒ 「4 笔业务」的具体划分没有书面来源。

## 3. 必须先问客户的 4 项（客户的决定，不能替她填）

| # | 问题 | 为什么影响结果 |
| :--- | :--- | :--- |
| 1 | 分析目的 (B22)：是否延续 P31（靶点 + 脱靶 + 载体）？有无新目标？ | 决定定制项范围 |
| 2 | Ctr_1–3 是 CRISPR-naive 亲本 RPE-1，还是带 Cas9 的处理克隆？ | 决定能否用 test − control 差集；P31 因对照不 naive 额外投入 6 h（§4） |
| 3 | Sample 1–7 各是哪个克隆/处理？gRNA、载体是否与 P31 相同？是否全为单克隆？ | 不同 gRNA/载体 ⇒ 脱靶与载体筛查要逐套做 |
| 4 | ~~浓度 `total 200ng/ul`：4 µg 还是总量 200 ng？~~ | **已关闭 2026-11-06**：建库已由 Hung 完成，问了不改变任何决定 ⇒ 见 §9 第 3 项 |

另：需新报价单号（表单沿用 QTE-26-07-30-001）；`Ctr_2 ` 名称带尾随空格。

### 3.1 测序深度：客户本批**未提任何要求**（2026-10-28 核查）

- 表单没有深度字段，B23 也是模板原文；Gmail 里 Yeji 的来信正文（不含引用部分）中，唯一和深度
  有关的是 08-28 那句 "would we need to do deeper sequencing?"（问的是 P31 靶点百分比）。
- 历史口径：P31 时她主动要 10×（2026-07-30 来信 "minimum of 10× coverage"，依据是介绍人
  postdoc 的意见；同时说这是 pilot，"deeper sequencing in the future maybe"）；Jianjun 当时也
  回过 "standard data analysis for coverage>=30x"、"higher-depth WGS (e.g., 30× or above)"。
- 我方 09-01 建议（见 §2）之后**客户未回复**。
- ⚠ 深度按规则是我方的决定，但测序价格按 X 计（Jianjun 报价公式：1–11 样本 45 + 11.8·X 美元/样本，
  ≥12 样本 45 + 10.8·X）⇒ 必须写进新报价，与 Jing 对齐。10 样本按公式：10× = $1,630，
  30× = $3,990，45× = $5,760（仅测序，不含定制分析）。

## 4. P31 的实际工时 —— 无计时记录

| 数字 | 性质 | 出处 |
| :--- | :---: | :--- |
| 定制部分 15–30 h | 开工前**预估**（2026-08-21），非实际 | P31 `docs/P31_research_plan_and_review_0821.md:292` |
| 因对照非 CRISPR-naive **额外 6 h** | 用户 2026-08-28 **事后认定**，非计时；已写入致客户邮件（"about six hours"），计费、折扣由 Jing 决定 | P31 `docs/with_client_info_20260825/P31_round2_plan_and_review_0825.md:173-180` |
| 最终开票小时数 | **未查到**（项目文档与 Gmail 均无发票记录，应在 Jing 处） | — |

**旁证：文件修改时间**（`scripts/ docs/ logs/ custom_research_report_2026082*/ client_materials/`，
`find -L ... -printf '%TY-%Tm-%Td %TH:%TM'` 按日取首末）：

| 日期 | 首次修改 | 末次修改 | 文件数 |
| :--- | :---: | :---: | :---: |
| 08-21 | 13:40 | 18:32 | 5 |
| 08-22 | 19:13 | 23:40 | 8 |
| 08-23 | 15:14 | 22:41 | 28 |
| 08-24 | 11:47 | 21:58 | 6 |
| 08-25 | 11:39 | 22:24 | 39 |
| 08-26 | 15:55 | 23:27 | 57 |
| 08-27 | 09:40 | 19:46 | 51 |
| 08-28 | 21:41 | 23:56 | 7 |
| 08-29 | 09:24 | 12:59 | 4 |
| 08-31 | 19:41 | 19:41 | 1 |
| 09-01 | 11:54 | 11:54 | 1 |

- 跨度 11 个日历日（08-21 → 09-01），集中工作约 9 天，每日首末窗口合计约 60 h。
- **⚠ 这是墙钟窗口，含机器计算、AI 执行和等客户回信的空档，只能当上限参考，不是人工工时。**
- 最重的 08-25 至 08-27：客户 08-25 才给 gRNA / 靶点 / 载体信息，导致整轮重做；P31 round2 文档 §6.3
  将其单列为「客户信息晚到导致的整轮重做」，**不并入**上述 6 h。

## 5. P46 工时估算（按 P31 范围推算，未经实际工时校准）

| 项目 | 类别 | 估算 (h) | 依据 |
| :--- | :---: | :---: | :--- |
| QC、比对、标准 germline / SV（Structural Variant）/ CNV（Copy Number Variation）、基础报告 | Regular | 0（不计费） | 默认管线 |
| 靶点 202 bp 缺失：逐样本接头读段、等位基因比例 | Custom | 3–4 | 复用 P31 `scripts/26_probe_client_reply_facts.py` |
| 脱靶：候选位点枚举（P31 为 1,587 个）× 10 样本逐位点证据，处理组 vs 3 对照 | Custom | 6–8 | 复用 P31 stage10；样本 2 → 10 |
| 载体整合筛查（Cas9-EGFP / mCherry-53BP1，先做同源屏蔽）及整合位点 | Custom | 5–8 | 复用 P31 stage11；须处理 53BP1 cDNA 假阴性陷阱 |
| 处理样本特有 SNV（Single-Nucleotide Variant）/ indel / SV / CNV（3 对照作基线） | Custom | 6–10 | P31 为 1 对 1，此为新统计设计 |
| 定制图表、报告撰写、客户沟通 | Custom | 6–8 | P31 往来 18 封邮件 |
| **合计 Custom** | — | **26–38** | 按 $120/h 约 **$3.1k–4.6k** |

**上浮条件**：
- 10 样本实为 4 笔独立业务（不同 gRNA / 目的）⇒ 按每笔分别估，可能明显高于 26–38 h；
- Sample 1–7 用了不同 gRNA 或载体 ⇒ 脱靶与载体筛查逐套做；
- 要求与 P31 两样本合并分析；
- 对照仍非 CRISPR-naive ⇒ P31 §6.1 所列 5 项额外工作照样发生（P31 round2 §6.5 建议计入报价）；
- 深度仍约 6× ⇒ 工作量未必增加，但大量结论只能写「无法判定」并需额外解释。

## 6. 下一步

- [ ] 客户回答 §3 四项后，重估 §5 并更新本文。
- [ ] 开工即**按 stage 记录起止时间**（P31 round2 §6.2 已立此要求），为后续报价提供实测依据。
- [ ] 若获 P31 最终开票小时数，补入 §4。
- [ ] 数据到达后：填 `scripts/0_sample_table.csv` 的 R1/R2 GiB、read pairs、深度，并在本文记录总量。

---

## 7. 客户问卷解读（2026-11-06，源文件 `client_materials/wgs_questionary.docx`）

读取方式：`python3 -I` 解包 `word/document.xml` 逐段提取（脚本在 scratchpad，非项目产物）。
问卷正文为 P31 当初那份问卷的同一模板，客户在每问下方作答。

### 7.1 客户给定的事实

| 项 | 客户答复（原文要点） | 与 P31 的关系 |
| :--- | :---: | :--- |
| 分析目的 | **全基因组脱靶分析**；靶点 ~200 bp 缺失；并"analyze any further mutations by the CRISPR system" | 与 P31 一致，加"其他突变" |
| 靶点 | PPP1R2，chr3；切点 **195,542,947 / 195,543,149**（间距 202 bp） | **与 P31 逐字相同** |
| gRNA1 | `gacttcctctatggtggcgt` + `cgg` | 同 P31 |
| gRNA2 | `gtttccgaggcagcagttg` + `cgg`（**仍写成 19 nt**） | 同 P31；19/20 nt 之问客户两轮都未正面答复 |
| 细胞系 | hTERT RPE-1 (**ATCC**)；全部 **单克隆** RPE-1_53BP1-mCherry；传代均 < 20 | 同 P31 亲本 |
| 外源序列 | **53BP1-FFR-mCherry（mCherry-BP1-2 pLPC-Puro）**，**逆转录病毒**包装递送，荧光分选后 3 cells/well 扩增 | **与 P31 同一载体**，图谱已在手（见 §7.3） |
| CRISPR 形式 | **dual caged gRNA**：caging 基团阻断切割，**光照解笼**，以此选择性地只在单个等位基因或两个等位基因上诱导活性 | P31 已是同一体系（"spatially controlled CRISPR"） |
| 对照 | Ctr_1–3 **含同一载体、同样经过 CRISPR 处理**，仅凝胶未见缺失；客户自述"not ideal as a control …… 但这是唯一的无单等位缺失单克隆" | **与 P31 同一处境** ⇒ P31 §6.5 的报价条件照样适用 |

### 7.2 样本编号：问卷里有两套编号，已对齐（偏移 3）

问卷中出现两种写法，**不矛盾**，相差一个常数 3：局部编号 `Sample 1–7` ↔ 全局编号 `4–10`
（全局 1–3 = Ctr_1–3）。两处交叉验证一致：

| 问卷出处 | 原文 | 解析 |
| :--- | :---: | :--- |
| Q1 / Q5 | 单等位 = Sample 1–4；双等位 = Sample 5–7 | 局部编号 |
| Q2 | "~50% 和 100% 效率见 sample 4–10，sample 1–3 无显著 CRISPR 活性" | 全局编号 |
| 末尾反问答复 | 单等位 = sample 4–7；双等位 = sample 8–10 | 全局编号 |

⇒ **结论**：`Ctr_1–3` 对照（无缺失）；`Sample 1–4` 单等位 Δ202；`Sample 5–7` 双等位 Δ202。
已写入 `scripts/0_sample_table.csv`（`client_global_no` 列即全局编号）。
**该对齐是我方推断（依据是两套编号在 7 个样本上的一一对应），需在致客户信中用一句话确认，
不得当成客户已确认的事实。**

### 7.3 "attached plasmap" 在本项目目录里不存在，但 P31 有同一构建体

- 本轮 `client_materials/` 实测只有 3 个文件（`Email.txt`、`Quote 10282602.pdf`、`wgs_questionary.docx`），
  **无质粒图谱附件**。
- P31 `client_materials/` 存有 **`mCherry-BP1-2 pLPC-Puro.dna`（57,736 B）**，构建体名与本轮问卷
  逐字相同 ⇒ **可直接复用**，不必再要一次。
- P31 另有 **`pMSCVpuro-Cas9-EGFP (2).dna`**。**本轮问卷未提 Cas9 载体**，只说细胞系"only contains
  53BP1-FFR-mCherry" —— 那么 Cas9 蛋白这一轮是怎么进去的（RNP 电转？瞬时质粒？还是这些克隆本身
  来自 P31 那个带整合 Cas9-EGFP 的系？）**没有答案**，这是 §9 第 1 项。

### 7.4 相对 P31 的三个实质变化（决定本轮分析设计）

| # | 变化 | 对分析的影响 |
| :--- | :--- | :--- |
| 1 | **深度 6.3× usable → 45× raw（约 27–30× usable）** | P31 大量"只能报有/无证据、不能报频率"的结论本轮可以定量；等位基因比例、杂合 indel、拷贝数都进入可判定区 |
| 2 | **样本 2 个 → 10 个（7 处理 + 3 对照）** | 从 1v1 并列升级为 **7 vs 3 组间比较**；"仅见于处理组、3 个对照全阴"是 P31 做不到的证据等级 |
| 3 | **新增双等位编辑组（Sample 5–7）** | 必须能区分「两个等位基因各自 Δ202」与「单等位 Δ202 + 另一等位更大范围丢失/LOH」——两者在凝胶上可能同样表现为单一条带 |

---

## 8. 分析步骤设计（Stage 0–5）

> 设计原则：**P31 round2 已跑通的判据与脚本能复用的全部复用**（同 gRNA、同参考、同载体 ⇒ 候选表与
> 屏蔽 mask 可直接继承），工时投在 P31 因深度/样本数做不到的部分。
> ⚠ 按 §2A：每个 Stage 在开跑前必须先把下表的「证据单位 / 保留与去重规则 / 判定指标 / 允许的结论」
> 落成脚本内的 decision table，并产出过滤前后计数的 audit 表；本节是该 decision table 的上层框架，
> 不替代逐 Stage 的实现级审查。

### 8.0 Stage 0 — 数据清点、比对、样本身份核查

| 项 | 做法 |
| :--- | :--- |
| 清点 | 复用 P31 `scripts/4_inventory_raw_data.sh`（RAW 改指本客户 Dropbox 目录）：文件数、大小、`gzip -t`、R1/R2 配对。**产出必须写成「截至 YYYY-MM-DD HH:MM 目录内为 N 个文件」**（Dropbox 渐进同步）。多 lane 时按 lane 汇总 |
| 不做 | **不跑 fastp/trimming** —— `Utilities/software_pitfalls/fastp.md` 记录的 1.3.3 数据驱动死锁；bwa-mem2 自带 soft-clip，接头不会伪造断点（假阴性方向，报告写明） |
| 比对 | `bwa-mem2 mem -t 36 -K 100000000 -Y` → `samtools sort -@ 8` → **GATK MarkDuplicates（非 Spark 版）** → index；**样本间串行**（36+8 = 44 线程 ≤ 56 上限） |
| 覆盖度 | **统一 mosdepth 默认模式**（含 mate-overlap 校正）`--by 1000`，全程不混用 `samtools depth`（配对重叠区可差 19%） |
| 身份核查 | 10 个样本同属一个单克隆亲本 ⇒ 全基因组公共位点基因型应高度一致。用多样本 mpileup 算两两一致率，**用于抓样本标签互换/污染**；一致率矩阵必须进报告（P31 为 2 样本版，本轮扩到 10×10） |
| 产物 | `analysis/stage0_<YYYYMMDD>/{bam,qc,cov,logs}`；回填 `scripts/0_sample_table.csv` 的 R1/R2 GiB、read pairs、raw/usable 深度 |

**资源与磁盘（按 45× raw 外推，须在首样本跑完后用实测速度回填）**：

| 量 | 估算 | 依据 |
| :--- | :---: | :--- |
| 单样本 raw | ~140 Gbp ≈ **45–55 GiB** gz FASTQ | 45 × 3.1 Gb 基因组；P31 10× 实测可按比例外推 |
| 10 样本 FASTQ | **~0.5 TiB**（落在 `/Work_bio`，实测可用 5.7 T） | `df -h /home/gao/Dropbox` |
| 单样本 markdup BAM | **~70 GiB** | P31 10× raw 实测 15.3 / 16.4 GB，线性外推 |
| 10 样本 BAM + 索引 | **~0.7 TiB**（落在 `/mnt/ex_8T_SSD`，实测可用 2.6 T） | 同上 |
| 临时 `sorted.bam` | 每样本额外 ~70 GiB（markdup 后即可删） | — |

⚠ **临时 `sorted.bam` 的删除需用户逐次同意**（CLAUDE.md「删除/覆盖前必须先获得用户同意」，
agent 无权自行判定"可再生中间产物"）。⇒ 计划里把它写成**显式请求项**：每样本 markdup + index 完成、
`samtools quickcheck` 通过后，列出路径与大小向用户请求一次批量授权；未获授权则按 10 × 140 GiB
预留磁盘（2.6 T 仍够，但会压缩其他项目空间）。

### 8.1 Stage 1 — 靶点解剖（PPP1R2 Δ202，10 样本）

- **复用**：P31 `scripts/with_client_info_20260825/2_stage9_ontarget.py` 的判据 ——
  断端读段（比对恰好终止于 cut1 / 起始于 cut2 且最大 soft-clip ≥ 20 bp）+ 接头认定（clip 序列与对侧
  下游参考一致度 ≥ 90%）。
- **本轮新增（深度带来的）**：等位基因比例 = 接头读段 /（接头 + 跨切点连续读段），30× 下杂合期望
  ~15 条接头读段 ⇒ **可以报比例，不只报有/无**。
- **本轮新增（双等位组）**：对 Sample 5–7 必须同时检验三种可能，不能只验证客户的预期：
  1. 两个等位基因各自 Δ202（期望：接头读段 ≈ 局部深度，跨切点连续读段 ≈ 0）；
  2. 单等位 Δ202 + 另一等位更大范围丢失（期望：缺失区两侧深度下降、断点在 202 bp 之外）；
  3. 单等位 Δ202 + 另一等位含切点处小 indel（凝胶分辨不出）。
  判据写在脚本里：**Δ202 区间内的归一化深度**（对该样本全基因组中位深度）+ 两切点各自的 indel 谱 +
  缺失区外 ±50 kb 的深度轮廓。
- **其他结局**：倒位（切点间片段反向回插）、串联重复、两切点之间的复杂重排 —— 用 clip 方向与
  mate 方向枚举，不预设只有缺失。
- **内建阳性对照 / 定标**：客户凝胶 + Sanger 已给出 10 个样本各自的预期基因型（`prior_evidence` 列）。
  **流程必须把 7 个处理样本判成缺失、3 个对照判成 WT**；任一不符 ⇒ 先怀疑流程或样本标签，
  不直接改写客户的实验结论。这是整条流水线唯一的外部真值，必须在报告里当校准点写出来。
- 产物：`analysis/stage1_ontarget_<YYYYMMDD>/`

### 8.2 Stage 2 — 全基因组脱靶（客户主诉求）

| 子步 | 做法 | 复用 |
| :--- | :--- | :---: |
| 2a 候选枚举 | 3 条 spacer 版本（gRNA1 20 nt、gRNA2 客户写的 19 nt、基因组补全的 20 nt `CGTTTCCGAGGCAGCAGTTG`）× 全基因组 ≤ 4 mismatch；**只用 spacer 搜，PAM 事后从参考实测**（不用 `seqkit locate -d`，避免 `-m`/`-d` 交互）；NGG 主表 + NAG 次表 | P31 `3_stage10_enumerate.sh`；**同 gRNA + 同参考 ⇒ 候选表可直接继承**，实测 `analysis/stage10_offtarget_20260825/tables/offtarget_candidates_with_evidence.tsv` 为 1,588 行（含表头）= **1,587 个候选位点**；继承前复核校验和并重跑一次定标断言（靶点必须被自己的流程检出），不重跑枚举 |
| 2b 逐位点证据 | 10 样本 × 每位点：indel、soft-clip 聚集、归一化深度、读段支持数 | P31 `5_stage10_inspect.py`，样本数 2 → 10 |
| 2c 成对切点脱靶 | 两条 gRNA 都够得着、间距 50–500 bp 的位点（P31 全基因组仅 6 个，全部是 PPP1R2 家族）；判据与 Stage 1 完全一致，**靶点本身是这一步的阳性对照** | P31 `13_stage10c_paired_cut_offtargets.py` |
| 2d 组间判定 | **7 处理 vs 3 对照**的样本级存在/缺失表 + Fisher 精确检验；事先写明 n=7/3 的检验力（最小可检出：7 阳性 0 阴性 ⇒ p = 1/120 ≈ 0.0083） | **本轮新增**（P31 只能 1v1 并列） |

**必须写进报告的边界**：
- 只支持 **mismatch**，不支持 DNA/RNA bulge（需 cas-offinder，本机未装；§4 规定未经确认不装依赖）⇒
  bulge 型脱靶属已知盲区；
- **对照不是 CRISPR-naive**（§7.1）⇒ "仅见于处理组"只能表述为**与编辑共分离**，不等于因果；
  同样地"处理组和对照组都有"不能一律判为背景 —— 对照也经过了 CRISPR 处理；
- 枚举只覆盖"序列近似"这一类脱靶，不覆盖脱靶后再经修复产生的远端事件 ⇒ 由 Stage 4 兜。

### 8.3 Stage 3 — 载体与递送物筛查

- **复用 P31 refs**（2026-11-06 实测确认内容）：
  `analysis/round2_refs_20260825/vectors.fa` —— **已含两条序列** `pLPC_mCherry_BP1_2` 与
  `pMSCV_Cas9_EGFP`（含 bwa-mem2 索引，可直接用）；
  `vector_human_mappable_mask.bed` —— **1 个屏蔽区间 `pLPC_mCherry_BP1_2:5550-7260`（1,711 bp）**，
  即载体内与人基因组同源的 53BP1 cDNA 段（150 bp 瓦片回贴 GRCh38，MAPQ ≥ 20 即屏蔽）；
  `vector_features.tsv` —— 载体注释（AmpR / ori / mCherry / EGFP 等），骨架阴性对照直接查这张表。
  ⚠ **这个 mask 是必需的，不是优化**：P31 第一版把 3 个"整合位点"全判在 **TP53BP1**，实为载体里
  53BP1 cDNA 与人基因组的同源假象（P31 round2 §5.2 坑 #8）。本轮载体相同 ⇒ 同一陷阱必然重现。
- 每样本：候选读段（unmapped / mate-unmapped / clip ≥ 20 bp）回贴载体 → 覆盖度 → 拷贝数；
  嵌合读段 → 整合接头位点（P31 坑 #6/#7：载体环状而参考线性，**clip 必须回贴载体而不是基因组**）。
- **内建对照**：P31 已实测两个整合接头位点（FADS1/FADS2 与 EMCN 附近），两样本断点完全相同 ⇒
  来自共同转导亲本。本轮 10 个克隆若同源于同一转导群体，**这两个位点应在 10 个样本中全部出现** ⇒
  当阳性对照用；出现**仅部分样本有**的新接头位点，才是克隆间差异。
- 骨架（AmpR / ori）出现即需解释（逆转录病毒转导不应整合骨架）——阴性对照。
- Cas9-EGFP 载体**已在** P31 的 `vectors.fa` 里 ⇒ 不论 §9 第 1 项如何答复，本轮都把它一起查：
  有信号则回答"Cas9 载体是否整合"，无信号则是一个干净的阴性结果。**但若本轮用的是另一个 Cas9
  构建体（而非 P31 那个 pMSCVpuro-Cas9-EGFP），则必须拿到新图谱**，否则该问题存在系统性盲区。
- ⚠ **计费边界**：载体筛查**不在**用户 10-29 实发邮件承诺的三项免费范围内
  （见 `docs/P46_email_depth_quote_1029.md`），但它是"CRISPR 造成的其他变化"的判读前提。
  是否计费由用户/Jing 决定，**分析侧照做并记录工时**。

### 8.4 Stage 4 — "CRISPR 造成的其他变化"（客户 Q1 末句）

这是本轮**设计上最需要先定口径**的一步（§2A）。10 个样本全部是**各自独立的单克隆扩增系**，
每个克隆在培养与扩增中会各自积累私有变异 ⇒ **逐个变异归因于 CRISPR 是不成立的**。
可成立的设计是把 3 个对照当作「同样处理、同样单克隆扩增、但无靶点缺失」的**本底率估计**：

| 要素 | 本 Stage 的口径 |
| :--- | :--- |
| 研究问题 | 编辑成功的克隆是否比未编辑的克隆携带**更多**的新发变异？若有，集中在哪类（SNV/indel/SV/CNV）、哪些位置？ |
| 证据单位 | **样本私有变异**（该样本有、其余 9 个样本均无），不是单个读段 |
| 保留与去重规则 | 按 read-name/fragment 去重计独立分子数；**但 primary / secondary / supplementary alignment 分开记录**（SV、融合、整合位点判读依赖 secondary/supplementary，不得因"不是独立分子"而丢弃 alignment context） |
| 判定指标 | ① 每样本私有变异计数的组间分布（7 vs 3，非参数检验 + 逐样本列出，不只报 p）；② 私有 SV/CNV 断点与 Stage 2 脱靶候选、Stage 3 整合位点的重合情况 |
| 允许的结论 | "处理组私有变异负荷高于/不高于对照组"；"某私有事件落在脱靶候选位点 ⇒ 有机制上的解释"。**不允许**把任一私有变异单独写成"由 CRISPR 造成" |
| 过滤审计 | 每条过滤给出过滤前后的 record / read-name / 候选数，保留 machine-readable excluded 表 |

- 工具：germline SNV/indel 用 10 样本联合 call（GATK HaplotypeCaller + GenotypeGVCFs，或
  DeepVariant 逐样本 + GLnexus，**开跑前查 `Utilities/software_pitfalls/INDEX.md`**）；
  SV 用 manta/delly；CNV 由 mosdepth 1 kb 窗 → 分段（复用 P31 `scripts/8a_cnv_segment.R`）。
- **缓存命名带样本数**（P41 教训）：联合 call 产物写成 `joint_n10.vcf.gz`，防止样本追加后静默误用旧缓存。
- 另：P31 第一轮的 16 个 SV 候选重解释 stage（P31 stage12）在本轮对应为「私有 SV 与脱靶/整合清单
  逐一比对」，已并入本 Stage，不单列。

### 8.5 Stage 5 — 交付

- 目录 `custom_research_report_<交付日>/`；报告结构沿用 P31
  `custom_research_report_20260828/Dual_Cut_Project_WGS_0828.md` 的 §1–§9（Objectives / Core
  Conclusions / Key Findings / Sample Information / Analysis Rationale / Methods / Results /
  Conclusions+Limitations / Deliverable Files）。
- 机械化自审：复用 P31 `11_check_report.py`（149 项）+ `12_key_metrics.py`，
  并按本轮新增内容扩项（10 样本的每样本数字必须对得上源文件；组间检验的 n 与 p 必须可复算）。
- 图件 PDF only；配色按 `scientific-figure-making` skill；交付前逐张人眼看图（P31 坑 #9–#13 全是
  只有人眼能发现的静默失效）。
- **交付物只准备、不发送**（CLAUDE.md「对外交付与署名」）。

### 8.6 执行纪律（开跑时必须照做）

| 项 | 要求 |
| :--- | :--- |
| 工时 | **按 stage 记起止时间**（P31 round2 §6.2 已立此要求，P31 本身没做到）⇒ 本轮是第一份可引用的实测工时 |
| tmux | 每个长 stage 一个 session，`tmux new-session -d -s <name> '... 2>&1 \| tee log.txt'`；用 `capture-pane` / tail 看，不交互式 attach |
| 总控 | 每轮**重新从 `scripts/0_sample_table.csv` 枚举**待处理样本，靠 `.done` 跳过；不在启动时 `ls` 一次定样本集合（Dropbox 渐进同步，P41 I-3 死锁教训） |
| 冒烟 | **每个消费脚本、每个样本类别都要冒烟**（P41 I-2/I-6）：本轮"样本类别"= 对照 / 单等位 / 双等位三类各取一个 |
| 早失败 | 启动后按 `watchdog` skill 做 Phase-1 早失败检测 |
| 归档 | 任何结果目录不 `rm -rf`；修正前旧版进 `superseded_for_comparison/<产物类型>/v<N>_<日期>_<缺陷>/` 并附 README |
| 开跑前 | 扫 `/home/gao/Utilities/software_pitfalls/INDEX.md`，覆盖本轮全部工具（bwa-mem2 / samtools / GATK / mosdepth / seqkit / manta / delly / DeepVariant），不只扫记得出过问题的 |

---

## 9. 待确认事项（**现行 3 项**：第 1、2、4 项；第 3 项已关闭）

| # | 问题 | 对象 | 为什么影响结果 | 若无答复的默认处理 |
| :--- | :--- | :---: | :--- | :--- |
| 1 | **本轮 Cas9 是怎么递送的？**（caged gRNA + Cas9 蛋白 RNP？瞬时质粒？还是这些克隆本身带整合的 Cas9-EGFP？）若涉及新载体，请给图谱 | 客户 | 决定 Stage 3 的 `vectors.fa` 组成与"骨架整合"阴性对照是否成立；若有 Cas9 载体而我方没纳入，整合筛查会有系统性盲区 | **必须问**。暂按"仅 53BP1-FFR-mCherry"建 refs，并把 P31 的 Cas9-EGFP 载体一并纳入作探查（有则报，无则为阴性对照），报告标注递送方式未确认 |
| 2 | **样本编号对齐确认**：Ctr_1–3 = 对照；Sample 1–4 = 单等位；Sample 5–7 = 双等位（§7.2） | 客户 | 判错组别 ⇒ Stage 2d/4 的组间比较整条失效 | 一句话确认即可；**不停等**，按 §7.2 的解析先跑（猜错只需改分组列、不必重跑比对） |
| 3 | ~~浓度 `total 200ng/ul` 是 4 µg 还是总量 200 ng？~~ | — | **已关闭（2026-11-06，用户裁定）**：Hung 已完成建库，测序即将出结果 ⇒ 该歧义不再影响任何可做的决定，**不问 Jing、不问客户** | 不追问。表单浓度字段按「客户书写原样」存档于 `scripts/0_sample_table.csv` 的 `conc_as_written` 列，`total_ug` 保持 `UNRESOLVED`；实际投入量以实验端记录为准。若测序 QC 出现文库复杂度偏低（重复率异常高、insert size 分布异常），再回头向实验端取实测定量**作为解释依据**，而不是作为决策输入 |
| 4 | **是否要与 P31 两样本合并分析**（12 样本一起看） | 用户（商务边界） | P31 深度仅 6.3×，与本轮 45× 混在一起会使"检出/未检出"不可比；且**合并分析不在 10-29 实发邮件承诺的免费三项内** | **不主动做**。本轮按 10 样本独立成篇；仅在报告里引用 P31 的阴性结论作背景 |

**已不必再问的两项**（§3 原第 1、2 项）：分析目的（Q1 已答）、对照性质（Q4 已答，确认为非 naive）。
**不再追问 gRNA2 的 19 vs 20 nt**：两轮未答，两版都做的成本很低（P31 已有脚本），
报告里并列呈现即可。

### 9.1 §5 工时估算的修订方向（不改数字，记录依据）

§5 的 26–38 h 是按 P31 范围 + 10 样本推的，当时未知的三件事现在已知，**全部指向上浮**：

| 已知事实 | 对工时的方向 |
| :--- | :--- |
| 深度 45× raw（而非 6×） | Stage 0 机器时间显著上升（人工时不一定），但**结论从"有/无"升级为定量 ⇒ 报告与校核工作量上升** |
| 建库已完成、测序即将出结果 | 投入量/文库参数已成既定事实，不再是待定变量 ⇒ **不产生额外沟通工时**；只有 QC 异常时才需回溯（§9 第 3 项） |
| 对照仍非 CRISPR-naive | P31 §6.1 的五项额外工作照样发生（P31 事后认定 6 h，样本数 2 → 10） |
| 新增双等位组 + 组间统计设计 | Stage 1 的三态判别与 Stage 4 的 7 vs 3 统计口径都是 P31 没有的新设计 |

⇒ **待 §9 第 1、3 项答复后一次性重估**，本轮按 stage 实测记时，使下一次报价有实测依据（§8.6）。
