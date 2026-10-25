# P45 核心运行记录与参数选择

- **创建：** 2026-10-25（修订记在文末更新记录）
- **范围：** 本文件记录执行参数、运行时间线和过程中的发现。范围决定与证据模型以 [`P45_research_plan_and_review_1025.md`](P45_research_plan_and_review_1025.md) 为准（SSOT），两边互相链接，不复制内容。
- **代码路径**均相对于 `scripts/20261025/`。`from_P44/`、`from_P42/` 是来源快照（2026-10-25 12:05，md5 见 `SNAPSHOT_MD5.txt`），本层的同名或同功能脚本可以直接 `diff`。

## 1. 影响结果的决定

| # | 做了什么 | 为什么 | 代码 | 独立复核 |
| :--- | :--- | :--- | :---: | :--- |
| 1 | Stage 1–8 不重跑，从 P44 第 2 轮（W ±50 kb）复制；Stage 9 从 P44 第 1 轮复制 | 计划 §2.2；用户要求可复用的直接拷贝 | `0_wait_and_import.sh` | `analysis/20261025_inherited/import_manifest.tsv` 中 md5_source == md5_dest |
| 2 | GRIDSS / SvABA 取 P42 的运行结果（GRIDSS 18 样本联合；SvABA 用 v5 DLD1 multisample，与 P42 `07_` 标准化所用版本一致） | 两者的 discovery window 与 W 完全相同（`discovery_window_0based.bed` = chr11:114553999–114675500） | `0b_import_external_callers.sh` | `cat analysis/20261025_inherited/P42/external_callers/discovery_window_0based.bed` |
| 3 | 8x 中 Manta/DELLY 用 P44 原函数重算，并断言与 P44 Stage 8 逐行一致 | 证明 4-caller 表的 Manta/DELLY 部分与 P44 完全相同 | `8x_caller_concordance_4callers.py`（文末 assert） | 用 P44 第 1 轮数据试跑：84,240 行一致 |
| 4 | GRIDSS 单断端（ALT 为 `.`）只记 one_end（anchor），不计 concordance | 单断端没有 partner | `svmatch.py` `pairs()` | — |
| 5 | caller 主口径为两端 ±100 bp + 方向一致；另输出 best_two_end_dist 与 ≤50 bp 标志 | 与 P41/P44 一致，同时能与 P42 的 ±50 bp 口径对上 | `svmatch.py` `classify()` | `concordance_long_4callers.tsv` 的列 `within_50bp` |
| 6 | claim ladder 按计划 §4.2 实现，结果出来之后不改判据 | 用户 Q2 批准 | `9x_structural_coherence.py` | `stage9x_coherence/coherence_tests.tsv`，逐样本列出判据 ①② |
| 7 | 客户可见的 junction ID 去掉内部前缀（P44J0332 → J0332）；列名 `P41_event` 改为 `event_in_20261009_report` | 客户文件中不出现内部项目编号（【P44←P42】B9） | `p45lib.py` `cid()`；`10_build_tables.py` | `grep -c P44J <交付目录>/tables/*.tsv` 应为 0 |
| 8 | G1：genome-wide SNV/indel，S01–S10 联合 calling，参数与 P44 的 chr11 calling 相同（`mpileup -q20 -Q20 --max-depth 1000`、`call -mv`、`norm -m -any`、QUAL ≥ 30）；只做主染色体 chr1–22、X、Y、M | 用户 Q3 = b；保持与 chr11 结果可比 | `G1_genome_wide_snv_indel.sh` | `summary.tsv` 中 chr11 记录数应等于 P44 chr11 VCF 的记录数 |
| 9 | G2：DELLY germline→somatic 流程。合并 S01–S05 的 discovery sites（`merge -y 0 -a 0.05 -v 5`，不预先过滤）；5 个 BAM 逐一 re-genotype；`delly filter -f somatic -t`（只打标签，不删记录） | 用户 Q4 = b；客户 A2 要 genome-wide 的 +Dox-only 事件 | `G2a_delly_dld1_regenotype.sh` | `G2_dld1_genome_wide/summary.tsv` |
| 10 | G2 的逐样本判据：present = 非参考基因型且 alt（DV+RV）≥ 3；absent = alt 0 且 ref（DR+RR）≥ 10；alt 1–2 = low_support；其余 = underpowered | absent 规则与 locus 表相同（AF 0.25 漏检 ≤ 6%） | `G2b_dld1_event_table.py` | — |

## 2. 运行时间线

| 时间（EDT） | 事件 |
| :--- | :--- |
| 2026-10-25 12:05 | 用户批准计划（全部按推荐）并授权全流程；登记 active project 与 claude_code autostart 条目 |
| 12:06 | G1 启动（tmux `p45_G1_snv`，632 个 5 Mb chunk，20 并发） |
| 12:07 | Stage 0 waiter 启动（tmux `p45_stage0_import`），等待 P44 第 2 轮 CHAIN_DONE |
| 12:09 | Stage 0b 完成：8 个文件，md5 不一致 0 |
| 12:10 | G2a 启动（tmux `p45_G2a_delly`）：合并后 12,953 个位点 |
| 12:23 | 分析链 waiter 启动（tmux `p45_chain`）：Stage 0 完成后依次跑 8x → 9x → 10 → 10c → 11 → 14 → 14b |
| 13:37 | P44 第 2 轮 CHAIN_DONE；Stage 0 开始复制（P44 脚本与 12:05 快照 md5 一致） |
| 14:36 | Stage 0 完成：282,382 个文件，6.8 G，md5 不一致 0 |
| 14:39–15:45 | 8x → 14b 全部完成（P45_CHAIN_DONE）；14 的 genome-wide grep 计数约 1 h |
| 15:21 | G1 完成：raw 10,104,360，QUAL ≥ 30 9,477,203；chr11 451,844，与 P44 的 chr11 call 完全相同 |
| 15:50 | 9x v1 → v2（新增 well_supported_junction_candidate），重跑 9x/10/10c/11；v1 已归档 |
| 15:51 | G2a 的 S05 单进程 genotyping 异常慢（约 1.2 MB/s，S04 约 19 MB/s）⇒ 并行启动按染色体拆分的 G2a2（先用 S04 验证与单进程结果逐条一致）；首次启动因脚本 bug 失败（见 §3），已归档后重启 |
| 15:58 | S05 单进程提速（约 55 MB/s），预计先于拆分版完成 |
| 16:09 | G2a 完成；停止备用拆分版并归档 |
| 16:10–16:22 | 发现并修正 SV 记录 END 丢失（v2 → v3 → v4），重跑 8x、10、10c、11、G2b |
| 16:25 | 第一次交付构建：硬检查拦下 `event_group` 列中的内部 ID 和自审中的坐标重名；修正后重建，SELF-AUDIT OK |
| 16:33 | Fig5 虚线下移，重建交付 + 自审 OK + `check_report.py` exit 0；git `f201558`（未 push） |

## 3. 过程中的发现

- **P42 对 C2 的 DELLY concordance 记录有误（内部；P42 交付包从未发出）。** P42 的 `candidate_manta_delly_concordance.tsv` 中，C2 在 S05 的 `delly_exact_two_end_record_count` 为 0。但 P41 的 S05 DELLY discovery BCF 中有 `BND00010059`（chr21:20,648,241 ↔ chr11:114,607,190，PRECISE，SR=3，FILTER=LowQual），与 C2 的两个端点（chr11:114,607,191；chr21:20,648,242）各相差 1 bp。复核命令：`bcftools view -H -i 'ID=="BND00010059"' /home/gao/projects_2026M7/41_Vedict_18_wgs/analysis/sv/delly/discover/S05_DLD1_c4_dox.delly.bcf`。P44 第 1 轮已正确记为 “DELLY (unfiltered)”。P45 用 P44 第 1 轮数据试跑时，同一事件得到 DELLY + GRIDSS + SvABA 三个 caller 的两端匹配，距离分别为 2、1、0 bp。
- 用 P44 第 1 轮数据试跑 9x：S04 的端粒 junction 与 S05 的 fold-back 判为 structurally coherent；S05 的 chr21 junction 距 S05 深度断崖 CI 约 1 kb，判据 ① 不满足，留在 limited；S03 的 10 bp inversion 不满足判据 ②（S03 没有 LOH）。说明判据有区分力。正式结果以第 2 轮数据为准。

- **claim ladder 空档**：计划 §4.2 没有为 Well_supported 定义层级。第 2 轮的两个 GM12878 共有翼侧 junction（J1365、J0533）因此被 v1 归为 limited。处理见计划更新记录 2026-10-25 16:00。
- **G2a2 脚本 bug**：`local c=$1 o="$O/parts/$c"` 中，`$c` 在 `local` 赋值之前就被展开，导致所有分片都叫 `.sites.bcf`，没有产生任何 genotype。已改为 `$1`。G1 与 G2a 的同类写法本来就用 `$1`，不受影响。失败的原始产物归档于 `superseded_for_comparison/analysis/v1_20261025_G2a2_split_empty_chrom_bug/`。
- **第 2 轮与第 1 轮的目标区结果**：4 个 ≥ 2 fragment 事件的坐标、before/after 计数和 AF 完全相同，只是 ID 重新编号（P44J0332 → J0922，P44J0115 → J0359，P44J0289 → J0786，P44J0389 → J1036）。目标区 junction 总数由 530 变为 527（Not supported 由 272 变为 269）。
- **SV 记录的 END 被丢失（P44 Stage 8 同样存在，影响 P44 尚未发出的交付）。** pysam 不在 `rec.info` 中暴露 `END`，P44 `pairs()` 用的 `info.get("END", pos)` 因此总是退回到 POS，所有 DEL/DUP/INV 记录被当成长度 0 的点事件。
  - 后果一，locus caller concordance：这类记录不可能与 junction 第二端匹配。修正后 256,392 行 Manta/DELLY 中有 148 行升级、0 行降级（`stage8x_concordance/manta_delly_diff_vs_P44_stage8.tsv`），全部在 ±50 kb 翼侧；目标区 caller 计数不变。
  - 后果二，G2：所有事件的长度为 0，Manta 匹配退化为单端匹配（+Dox-specific 中 "with Manta" 的计数由错误的 6 变为正确的 7）。
  - 修正经过两步：v3 用 `rec.stop`，但 DELLY genotyped BCF 的 rlen 写为 1（例：DEL00003167 实际 END = 157,172,208，而 `rec.stop` = POS）；v4 改为从 INFO 文本解析 `END=`。v2、v3 产物分别归档于 `superseded_for_comparison/analysis/v2_20261025_svmatch_END_ignored/` 与 `v3_20261025_svmatch_END_from_rec_stop/`。
  - 独立复核：`bcftools view -H -i 'ID=="DEL00003167"' analysis/20261025_p45/G2_dld1_genome_wide/dld1_S01_S05.somatic_tagged.bcf | cut -f8 | tr ';' '\n' | grep ^END`
- **G2 结果**：12,953 个事件中 13 个为 +Dox-specific（S04 8、S05 5），两个克隆共有的为 0；7 个也被 Manta 调出；5 个通过 DELLY somatic filter。S05 单进程 genotyping 最终于 16:09 完成（中段异常慢的原因未查明）；备用的拆分版未使用，已归档。

## 更新记录

- 2026-10-25 — 创建。
