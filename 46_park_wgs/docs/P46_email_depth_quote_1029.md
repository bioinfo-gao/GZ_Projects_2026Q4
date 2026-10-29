# P46 — 致 Yeji Park：建议 45–50× raw WGS + 附 Quote 10282602（2026-10-29 草稿）

- **状态（2026-10-29 回读 Gmail Sent 核实）**：**用户已修改并发出**，2026-10-29 14:52:10 UTC，
  Sent message id `1a188259fa838bdc`，线程 `19fdd85e1110e93e`（回复我方 09-01 那封 `1a05e4746f78e525`）。
  附件 `Park_Quote 10282602.pdf`（用户挂的，文件名与本机存档不同）。**实发正文见文末「实发版」一节**；
  下面「草稿正文」是我方 Gmail 草稿 `r-1512917762955993518` 的原文，留作对照。
- 本机 quote 存档 `client_materials/Quote 10282602.pdf`，480,586 B，sha256[:16] `4d89c138b985ff62`
  （来源：用户 2026-10-29 14:38 UTC 发给自己的 "Park quote" 邮件）。
- **Quote 内容**（pdftotext 读取）：Quote 10282602，Created 10.28.2026，WGS 45x × 10，unit $576，Total $5,760。
  与 Jianjun 报价公式 45 + 11.8 × 45 = 576 一致。
- **To**：Yeji.Park@childrens.harvard.edu；**Cc**：j.shi@athenomics.com、kchung@athenomics.com、jianjun@athenomics.com
- **依据**：09-01 实发版（30x usable / 45–50x raw）；P31 实测 raw → usable 比 6.33/10.6 = 0.60、6.14/9.8 = 0.63
  ⇒ 45–50× raw ≈ 27–31× usable，信中写 "about 27x to 30x"。时间线依据 Kuo-Chan 2026-10-27 来信
  （本周建库、下周测序）。

## 草稿正文（我方原稿，未发出）

```
Subject: Re: Inquiry About Whole Genome Sequencing Services [EXTERNAL]

Hi Yeji,

Thank you for sending the 10 DNA samples for the next round. Please find attached our quote (Quote 10282602) for whole-genome sequencing (WGS) of all 10 samples at 45x raw output.

Following my email of September 1, we recommend 45x to 50x raw output per sample. Raw output means the total sequence produced divided by the genome size. In the first round, duplicate reads and the overlap between the two reads of each pair reduced about 10x raw output to about 6.3x usable depth, so 45x to 50x raw output gives about 27x to 30x usable depth. The quote is priced at 45x, the lower end of this range.

The reason is the one in my September 1 email: off-target detection, and any comparison between samples - structural variants, copy number, small insertions and deletions - depend on depth. At 6.3x, a change carried on one of the two alleles is covered by about three reads, which can fall to zero by chance. At about 30x it is covered by about fifteen reads, so a negative result becomes a real negative and a difference between samples becomes meaningful.

The libraries are being prepared this week and sequencing is scheduled for next week, so this is the point at which the depth is set. Please let us know if you have any questions about the quote.

Best regards,

Zhen Gao, PhD
Principal Bioinformatics Scientist
Athenomics
```

## 实发版（2026-10-29 14:52:10 UTC，Gmail Sent `1a188259fa838bdc`）

相对草稿的改动（用户所改）：
1. **第 2、3 段之间新增一段 —— 对客户的价格承诺**（原文照录）：

   > With sequencing at this depth, the data analysis for the questions you set out in the first
   > round - the deletion at the PPP1R2 target site, the genome-wide off-target search, and other
   > changes caused by the CRISPR system - will be provided for all 10 samples at no additional charge.

2. 署名下加电话 617-308-8078。其余各段与草稿逐字一致。

⚠ **免费范围以上面这句为准，后续报价/计费必须据此划界**：
- **免费**：10 个样本的 ① PPP1R2 靶点缺失；② 全基因组脱靶搜索；③ CRISPR 系统造成的其他变化。
- **未承诺免费**（信中没有提到）：载体整合筛查（P31 时我方加做的）、P46 与 P31 数据合并分析、
  客户之后新提的问题。
- 实发版**没有**写「超出部分另行报价」这句边界句（建议稿里有，用户未采用）⇒ 若客户之后把新需求
  归入「CRISPR 造成的其他变化」，只能逐项按上面三条原文判定。
