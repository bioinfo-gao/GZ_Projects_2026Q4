# P44 核心运行记录与参数选择

- **创建：** 2026-10-25（修订记在文末更新记录）
- **范围：** 本文件记录**执行参数**与运行时间线；客户要求、范围决定与证据模型在 [`P44_research_plan_and_review_1024.md`](P44_research_plan_and_review_1024.md)（SSOT，两边互相链接不复制）。
- 代码路径均相对于 `scripts/`；行号于 2026-10-25 01:40 grep 回对。复核命令尽量**独立于分析代码**。
- 正式输出：`analysis/20261024_read_rescue/`（stage2–11）；6 样本试跑在同目录 `*_test`（含 `stage5_assembly_test_v1_spadesfail`、`stage6_genotype_test_v1_nogroups` 两个试跑历史版），**不进交付**。

## 1. 影响结果的决定

| # | 做了什么 | 为什么 | 代码 | 独立复核 |
| :--- | :--- | :--- | :---: | :--- |
| 1 | 抽取窗口 W = target ±5 kb = chr11:114,599,000–114,630,500 | > 最大 insert；覆盖 114,601.5–603.5 kb 低复杂度段与深度边界 | `0_common.sh:20`、`p44lib.py:18` | — |
| 2 | 全 BAM 流式扫描，收集 mate / SA / XA 落在 W 的 record（regex 为 W 的超集，Stage 3 精确复核） | 原 BAM 无 secondary record，多重比对只在 XA tag | `2_genome_scan_rescue.sh:13-14` | `samtools view -c -f 256 <任一 P41 BAM> chr11:114599000-114630500` → 0 |
| 3 | 断言 W 内无 secondary、无 hard clip、无缺 SEQ（18/18 通过） | `SEQ` 完整性审计（§2A） | `1_window_inventory.py`（assert 段） | `analysis/.../stage1_inventory/class_summary.tsv` |
| 4 | soft clip ≥ 20 bp 才作断点证据 | 与 P41 口径一致 | `p44lib.py:20` | — |
| 5 | Adapter read-through：3′ segment 起点 ≥ adapter 起点 − 5 的 split 不计；adapter 在 clip 前 10 bp 内的 3′ clip 不计 | insert 中位 ~160 bp < 150 bp read；poly-G 尾比对到 chr2:32,916,3xx | `4_evidence_clusters.py:65`、`:121`；`p44lib.py:223` | `samtools faidx <REF> chr2:32916330-32916480` 为纯 G；`samtools view <S04 BAM> chr2:32916340-32916360` 可见 insert+`AGATCGGAAGAGC`+index+polyG |
| 6 | Low-complexity（entropy < 1.2 bits）segment / clip / discordant mate 不计 | poly-G 类 artefact | `p44lib.py:220`；`4_evidence_clusters.py:68,139,168` | — |
| 7 | Secondary placement 仅在 AS ≥ primary AS − 10 时作证据 | secondary 的科学含义是"竞争 placement"；`-a` 下大量远低于 primary 的局部 hit（S04 W 内 secondary 66,300 条） | `4_evidence_clusters.py:35` | `stage4_clusters/filter_impact.tsv` 的 `realn_secondary_(not_)competitive` |
| 8 | Discordant 判定：非 FR 或 \|TLEN\| > max(500, p99.9) | 按样本 insert 分布 | `4_evidence_clusters.py:201` | `stage4_clusters/maxins.tsv`（18 样本均 500） |
| 9 | Cluster 容差 15 bp（与 P41 同） | — | `4_evidence_clusters.py:33` | — |
| 10 | 局部组装只对 ≥ 2 非 dup fragment 的 cluster；SPAdes `--sc --only-assembler --cov-cutoff off -k 21,33,55`，read 先在 adapter 截断（≥ 30 bp） | 单 fragment 组装不添独立证据；微小 read 集触发 k-mer coverage model 报错（exit 64） | `5_local_assembly.py:56,116` | `stage5_assembly/assembly_log.tsv` |
| 11 | 竞争性基因分型：模板 400+400 bp；read-end 跨 junction 每侧 ≥ 20 bp，且 AS_alt ≥ AS_genome + 10 且 ≥ AS_ref + 10 | "after rescue" 的定义；同一 repeat 拷贝能同样解释的 read 不计 | `6_junction_genotyping.py:40,229` | 抽查 `stage6_genotype/read_support.tsv.gz` 中 P44J0332 的 7 行（AS_template vs AS_genome） |
| 12 | Event grouping：共享 ≥ 50% alt read-end 的 junction 合并为一个 event，计数取并集 | partner 位于多拷贝序列（端粒重复）时避免一事件多计 | `6_junction_genotyping.py:310` | `stage6_genotype/junctions.tsv` 列 `alternative_representations` |
| 13 | absent = 0 alt 且 ≥ 10 ref-spanning fragment，否则 underpowered | P41 口径 | `6_junction_genotyping.py:372` | — |
| 14 | Flank 唯一性：junction 两侧 100 bp `bwa-mem2 mem -a`，AS ≥ best − 10 的 placement 数 | §4 判据 | `7_breakpoint_annotation.py:40` | `stage7_annotation/flank_alignments.tsv` |
| 15 | 分级（§4 + Limited） | 见计划 §4 与 2026-10-25 更新记录 | `7_breakpoint_annotation.py:228-238` | — |
| 16 | Caller 匹配：两端 ±100 bp 且方向一致 | 与 P41 一致 | `8_caller_concordance.py:32` | `bcftools view -H <S05 delly bcf> chr11:114607000-114607400` 可见 chr21 BND |
| 17 | chr11 SNV：`bcftools mpileup -q 20 -Q 20` + `call -mv`，QUAL ≥ 30；BAF 位点 = 基线全部样本 AF 0.25–0.75 且 DP ≥ 10 | Q3：只做 chr11 | `9b_chr11_snv_calling.sh:15`、`9b_chr11_baf.py:40` | — |
| 18 | 深度边界：200 bp bin，以本样本 chr11:114.40–114.58 Mb 中位数归一，再除以 S01–S03 均值；two-level step + 1,000 block bootstrap（seed 44） | 给 junction 一个独立位置约束 | `9c_depth_boundary.py:24-25` | `samtools depth -a -Q 20 -r chr11:114605000-114606000 <S04 BAM>` 目视 |

## 2. 运行时间线

| 时间（EDT） | 事件 |
| :--- | :--- |
| 2026-10-24 23:19 | 新建 conda env `p44_assembly_20261024`（SPAdes 4.3.0），EXIT=0 |
| 23:22 | Stage 2 扫描 18 BAM（6 并发，每样 ~25–32 min）启动，tmux `p44_scan` |
| 23:36 | chr11 SNV joint calling（S01–S10），tmux `p44_snv`；~00:14 完成（raw 478,856 / QUAL30 451,844 records） |
| 23:41 | 链式 Stage 3→8 启动，tmux `p44_chain` |
| 23:50–00:50 | 6 样本试跑暴露并修正：secondary SEQ='*'；adapter/poly-G artefact；secondary 竞争性；SPAdes exit 64；event grouping（见计划更新记录） |
| 00:50 | Stage 2/3 全部 18 样本完成；Stage 3b 启动 |
| 01:12 | Stage 4–8 完成（CHAIN_DONE, EXIT=0） |
| 01:20–01:40 | Stage 10 表、Stage 11 图、报告源、交付包构建、self-audit |

## 3. 事故与误报

- `rm -rf` 被 hook 拦：目标是 scratchpad 里的 SPAdes 测试目录，但命令行含 `analysis/` 路径；改用新目录名，未删除任何文件。
- 9c 首跑 bootstrap 索引长度差 1（block 被末端截断）→ 改为只抽完整 5-bin block。
- 8 首跑 pysam 对 header 未声明的 INFO key 抛错 → header-safe getter。
- 交付报告首版 self-audit 抓到 2 处：S04 左侧深度比 0.985 写成 0.99（应为 0.98）；"Twelve" 与表格数字形式不一致。已改源文件后重建。

## 4. 实测结果（回填）

- 目标区 junction 530 个：Well 0 / Limited 3 / Ambiguous 1 / Weak 254 / Not supported 272（其中 6 个有 before ≥ 1）。
- 关键 event：P44J0332（S04，chr11:114,605,216–218 → (TTAGGG)n，0→4 fragments）；P44J0115（S05，376 bp fold-back）；P44J0289（S05，= P41 EV0623 = P42 C2）；P44J0389（S03，10 bp inversion）。
- 深度边界：S04 114,605,200（CI 605,000–605,400）；S05 114,605,800（CI 605,600–606,000）。
- Distal chr11q LOH：S04 0.978、S05 0.980；同一 homolog 保留 99.78%（8,004 位点）。
- 与 P42（内部，从未交付）对照：C2 = P44J0289；C3 = P44 Weak（S05 1 fragment）；C1（114,564,286）、C4（114,598,286）在 W 之外（P42 窗口 ±50 kb），属范围差异；P42 未报告 P44J0332 与 P44J0115。

## 更新记录

- 2026-10-25 — 创建。
