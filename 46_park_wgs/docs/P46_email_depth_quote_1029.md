# P46 — 致 Yeji Park：建议 45–50× raw WGS + 附 Quote 10282602（2026-10-29 草稿）

- **状态**：Gmail 草稿（draft id `r-1512917762955993518`，message id `1a1881cb00fea764`），**未发送**。
  回复我方 2026-09-01 那封（message id `1a05e4746f78e525`，线程 `19fdd85e1110e93e`）。
- **⚠ 附件未挂上**：Gmail 连接器无法从本机挂 PDF；用户发出前需手动附 `Quote 10282602.pdf`
  （本机存档 `client_materials/Quote 10282602.pdf`，480,586 B，sha256[:16] `4d89c138b985ff62`；
  来源：用户 2026-10-29 14:38 UTC 发给自己的 "Park quote" 邮件）。
- **Quote 内容**（pdftotext 读取）：Quote 10282602，Created 10.28.2026，WGS 45x × 10，unit $576，Total $5,760。
  与 Jianjun 报价公式 45 + 11.8 × 45 = 576 一致。
- **To**：Yeji.Park@childrens.harvard.edu；**Cc**：j.shi@athenomics.com、kchung@athenomics.com、jianjun@athenomics.com
- **依据**：09-01 实发版（30x usable / 45–50x raw）；P31 实测 raw → usable 比 6.33/10.6 = 0.60、6.14/9.8 = 0.63
  ⇒ 45–50× raw ≈ 27–31× usable，信中写 "about 27x to 30x"。时间线依据 Kuo-Chan 2026-10-27 来信
  （本周建库、下周测序）。

## 正文（与草稿逐字一致）

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
