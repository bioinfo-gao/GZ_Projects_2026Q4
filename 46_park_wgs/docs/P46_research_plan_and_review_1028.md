# P46 — Yeji Park 第二批 WGS（10 样本）：表单解读、待问清单与工时估算

- **Plan Date（创建，不可更改）**：2026-10-28
- **项目目录**：`/home/gao/projects_2026Q4/46_park_wgs/`
- **客户**：Yeji Park，HHMI / Boston Children's Hospital（TaekJip Ha lab）
- **前序项目**：P31 `/home/gao/projects_2026M7/31_Yeji_Park_RPE1_CRISPR_dual_cut_wgs/`
  （hTERT RPE-1，PPP1R2 双切 202 bp 单等位缺失 + 全基因组脱靶 + 载体整合筛查）
- **客户表单**：`Sample_Registration_Form_Yeji Park.xlsx`（项目根目录）
- **样本表（SSOT）**：`scripts/0_sample_table.csv`（本文不重复样本清单）
- **当前状态（2026-10-28）**：只有表单，**无测序数据**；分析范围未定，**未开跑**；未建邮件草稿。

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
| 浓度 (D42–D51) | `total 200ng/ul`，体积 20 µL | **有歧义**，见 §3 第 4 项 |

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
| 4 | 浓度 `total 200ng/ul`：200 ng/µL × 20 µL = 4 µg，还是总量 200 ng（≈10 ng/µL）？ | 后者低于表单要求的 2–10 µg，影响建库方式 ⇒ 转 Jing 确认 |

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
