# 归档：1127 报告在补入 Sample_5 第二条等位分子细节之前的版本

- **归档日期**：2026-11-27
- **归档文件**：`Dual_Cut_Project_WGS_1127.md`（与当日 11:52 交付目录内的版本逐字节相同）
- **这一版的缺陷**：
  1. §7.1 表与 §8.1 表把 Sample_5 第二条等位只写作 "a small indel at cut 1"，
     **未给大小与坐标**，而该 indel 的实测值（8 bp，`chr3:195,542,947–195,542,954`）
     就在同一交付包的 `on_target/indels_near_cutsites.tsv` 里；
  2. §8.4「Suggested next steps」写着
     `What is the second allele in Sample_5? → Sanger-sequence across cut 1 to read the exact indel`
     —— **建议客户去做一个已交付数据已经回答的实验**。
- **怎么改的**：按用户 2026-11-27 决定，补入坐标 + 移码判定 + 注释层面的 45 aa 截短预测，
  并把 §8.4 那条改为「Sanger 作为独立确认」而非「用来读出是什么」。
- **并排对照**：`diff` 本目录的 .md 与
  `../../../custom_research_report_20261127/Dual_Cut_Project_WGS_1127.md`。
- **证据与推导过程**：见
  `../../../Biology/PPP1R2_locus_and_RPE1_line_identity_1127.md` §5.3（含 CIGAR 复算命令）。
