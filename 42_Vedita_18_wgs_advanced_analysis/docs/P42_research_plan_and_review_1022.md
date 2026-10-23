# P42 全队列重复区 SV / CNV 证据重建 — 研究计划与审阅记录

- **首次建立：** 2026-10-22
- **本次实质修订：** 2026-10-23
- **状态：** 用户已完成旧版审阅并授权正式重跑；本文件是当前唯一有效的计划。
- **范围：** 全部 18 个 WGS 样本（S01–S18）。DLD1 的 S04、S05（+Dox）与 S01（parental）、S02/S03（control clones）构成主要因果比较；其余 13 个样本用于队列背景、技术/重复区可比性和细胞系特异性解释，不能被错误地当作 DLD1 的生物学对照。
- **核心位点：** `chr11:114,604,000-114,625,500`（GRCh38/hg38；1-based closed interval）。发现窗口预注册为其两侧各扩展 50 kb：`chr11:114,554,000-114,675,500`。
- **规范样本表（SSOT）：** `scripts/0_sample_table.csv`。
- **历史版本：** 旧的 S01–S10、P41-candidate-gated、仅重比对工作流及实际产物原样复制在 `superseded_for_comparison/P42_candidate_gated_read_rescue/v1_20261022_incomplete_window_evidence/`；它们保留为可比较原始记录，**不是**本轮主线结果。

## 1. 科学问题与结论边界

已观察到 S04/S05 的覆盖度低于 DLD1 非 Dox 样本，客户的 FISH 亦提示存在染色体片段异常。这是需要严格检验的线索，而不是 SV 已定位或 FISH 已由 WGS 独立证实的结论。本轮逐项回答：

1. S04/S05 是否在 chr11 靶位点或全基因组存在与缺失、重排或易位一致的 copy-number / breakpoint 证据？
2. 若有 junction，两个端点坐标、方向、microhomology 与 insertion sequence 是什么；其定位是否因重复序列而存在同等合理的替代位置？
3. 该事件在 S01/S02/S03 中是否缺失或显著低频，因而可称为 Dox-clone-specific；还是本来就存在于 parental/control？
4. 覆盖度降低是一段连续 CNV、全样本深度/文库差异、可比对性下降，还是多因素组合？

结论随证据等级而定：只有同时通过分子去重、断点两端证据、替代定位审计和对照比较的事件才称高置信。单端 soft clip、MAPQ0 或 caller 标记只能称候选/模糊。FISH 与 WGS 是独立证据类型，不能相互替代。

## 2. 输入、版本与不变性

- 输入 BAM：P41 mark-duplicate BAM/BAI，`/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs/analysis/bam/`；仅只读查询，绝不全基因组重比对、编辑或复制 BAM。
- 输入 SV：P41 Manta 与 DELLY 的 unfiltered/filtered 调用及索引；本轮重新解析记录和原始 read evidence，不把旧汇总表当真值。
- Reference：P41 的 GATK GRCh38 analysis-set；local-assembly contig 和 rescue reads 保持同一 reference，并保留所有替代 placement。
- 输入量：18 个 BAM、约 846 GiB（操作性估计）。新输出仅写入 `analysis/20261023_fullcohort_repeat_resolve/`，日志仅写入 `logs/20261023_fullcohort_repeat_resolve/`。重跑或参数改变创建新阶段/时间戳子目录，绝不覆盖旧结果。

## 3. 已发现的旧版问题与修正

旧版以 P41 的 179 个候选和已记录 QNAME 为入口，并使用 `read.query_name not in wanted or read.is_secondary`。这不是靶窗口的完整证据模型：遗漏了不在 179 个候选中的 read，也在比较替代位置前直接排除了 FLAG `0x100` secondary alignment。对于重复区“真实位置”的问题，这造成确认偏差。

本轮不再用候选表限制 read universe。discovery window 内全部 **primary、secondary（0x100）、supplementary（0x800）** records 都保留并分层；MAPQ0/multi-mapping、soft-clipped records，以及“mate 在窗口而本身 unmapped”的 records 全部纳入。`SEQ='*'` 仅表示该 BAM record 未存储 sequence；它被单列为“不能 sequence rescue”，绝不被误称无证据或 unmapped。

分析单位明确分两层：

- **alignment-record 层：** 保留 mapping hypotheses，判断重复区定位歧义；record 数不能作为独立支持数。
- **molecule 层：** 以 read name + mate/read-end 归并，记录 primary/secondary/supplementary 所有成员；同一分子至多给一个独立支持。duplicate flag 作为属性保留，在主计数与敏感性分析分别报告。

每一个 filter 输出 filter-impact audit：输入/输出 record、QNAME、molecule、sequence-available 数与原因。P41 179 行仅作为“新全窗口是否重新发现/支持”的注释，绝不是发现闸门。

## 4. 预注册分析设计

### 4.1 全窗口 read rescue 与局部重建

1. 每样本枚举 discovery window 全部 alignment records，并从窗口内 mapped record 追取其 unmapped mate。写出 record 清单、molecule ledger 和 FLAG/MAPQ/sequence-availability 审计。
2. 每个分子保留所有 sequence-bearing records；sequence 缺失及 hard/soft clipping 单列。MAPQ0、secondary、supplementary 不因类别丢弃。
3. 以 soft clips、split/SA、discordant pairs、跨断点 contig 和 caller 端点聚类生成**新候选**；先检查每个 read 的全 mapping 集，再评估对 chr11 的支持。
4. 每样本、每候选独立 local assembly（不同样本 raw reads 不混合）。contig 允许多重 placement 地重新比对 GRCh38；只有双端锚定、方向一致且替代 placement 已审计时才报告 base-pair junction。microhomology/insertion 从 contig/read CIGAR 与局部序列比对导出，保存序列来源。
5. 每候选分别给出 standard high-confidence 与 inclusive-rescue 计数：前者为主规则下独立 molecule；后者显示 MAPQ0/secondary/supplementary/mate rescue 贡献。两者绝不混写为“read count”。

### 4.2 SV caller 与全基因组结构背景

重新读取每样本 Manta/DELLY raw 与 filtered records，按 breakpoint-end、方向和 SV type 标准化，生成 `Manta_only`、`DELLY_only`、`both_concordant`、`both_nonconcordant`、`neither_caller_but_read_assembly_candidate`。caller concordance 提供独立算法一致性，不是高置信 junction 的充分条件；重复区无 caller 记录也不自动否定 read-level evidence。

全队列做两种分辨率 coverage/CNV 审计：全基因组 100-kb bins 检查样本深度和大片段 CN；chr11 10-kb、靶位点 1-kb bins 定位连续 CN boundary。standard CN 使用 primary、non-duplicate、非 QC-fail、非 supplementary、MAPQ≥20 records；另给 MAPQ≥0 inclusive track。每轨明确 flag、MAPQ、normalisation 与 mappability 限制，禁止把重复区 mapping-depth 当真实 CN。

### 4.3 DLD1 clone-specific 规则

主比较为 S04/S05 对 S01/S02/S03；S06–S18 仅作细胞系背景。候选合并条件为同 SV type、同方向、两个端点各在 50 bp 内；base-pair assembly event 以精确 sequence 优先。每候选列出五个 DLD1 样本的独立 molecule、总可评估 molecule/局部 depth、standard/inclusive AF、record-class 构成及 caller support。

- **Dox-shared high confidence：** S04 和 S05 各 ≥3 independent molecules，至少一个样本有双端/assembly junction；standard AF ≥5%；S01/S02/S03 同 junction 各为 0 或低于预计算检测下限。
- **Dox-private high confidence：** 仅 S04 或 S05 满足上述 molecule/AF/双端条件，全部 DLD1 对照不满足；只称 clone-private，不泛称 Dox-induced。
- **candidate / ambiguous：** 2 molecules、单端证据、MAPQ0/secondary 主导或多个等价 placement；保留、图示、计数，但不宣称定位明确或 clone-specific。
- **negative at available sensitivity：** 对没有观察的样本，给有效覆盖度和检测下限；不写“事件不存在”。

`≥3 molecule / ≥5%` 是初始阈值，原始计数始终并列报告。若 5% 小于 3 molecule，则 3 molecule 仍为下限；局部有效深度不足时降级为不可判定。

### 4.4 SNV/indel 与交付

对 DLD1 五样本核验同一 reference、同一过滤定义下的 SNV/indel VCF；若已有 P41 callset 的 provenance/完整性不合格，才在新目录补充 calling。全 18 样本 BAM/BAI 和 Manta/DELLY VCF 是既有交付源；本轮生成 manifest、校验和与 annotated per-event TSV，不复制 846-GiB BAM 到 P42 工作区。

## 5. 图片交付（PNG、PDF、来源 TSV 与脚本齐全）

至少制作下列 12 组图，且所有图片均可从保存的 TSV/脚本重建：

1. 全 18 样本 genome-wide normalized coverage heatmap/overview；
2. S04/S05 相对 S01–S03 的 chromosome-level coverage 差异；
3. chr11 10-kb CN tracks（standard / inclusive 并列）；
4. 靶窗口 1-kb depth、MAPQ、有效 molecule tracks；
5. 全样本 FLAG / MAPQ / sequence-availability audit；
6. read rescue 前后 independent-molecule 计数与类别构成；
7. Manta–DELLY–assembly/read evidence caller-concordance 图；
8. DLD1 五样本候选 junction AF/support heatmap；
9. 每个重点 junction 的端点、方向、microhomology 示意；
10. 每重点事件 read/contig alignment evidence panels（含替代 placement/歧义注记）；
11. local-assembly contig 双端锚定与替代位置比较；
12. FISH 观察、coverage、caller、junction evidence 的支持/未支持边界总结图。

图用色盲友好、白底可读的固定样本配色；不使用 pie chart。若没有 high-confidence junction，图 9–11 仍展示最强候选及不能定位的原因。

## 6. 执行步骤与质量闸门

1. `scripts/20261023/0_common.sh`：路径、reference、区域、全样本断言。
2. `1_inventory_full_window.py`：完整 record/molecule/mate rescue inventory。闸门：18 个完成标志；secondary/supplementary/MAPQ0/SEQ=* 都有显式计数。
3. `2_make_coverage_bins.py` 与 `3_coverage_cnv.sh`：两套 coverage。闸门：每 bin 可追到样本、flag/MAPQ 口径，样本总深度与 BAM 规模相符。
4. `4_parse_callers.py`：Manta/DELLY raw/filtered 标准化和匹配。闸门：明确每份 caller 文件是否存在、可读及 record 数；空文件绝不伪装阴性。
5. `5_cluster_junction_evidence.py`：从全窗口、非 P41 candidate list 聚类断点并写 molecule 级支持。
6. `6_local_assembly_and_align.sh`：按样本/候选 assembly、全 placement 比对、junction 判定。
7. `7_dld1_comparison.py`：DLD1 阈值、AF、检测下限与敏感性分析。
8. `8_make_figures.py` / `9_write_report.py`：图片、报告与主张→来源对照。

超过数分钟的步骤仅在 tmux 运行，单一长任务最多 28 physical cores / 56 threads，并启动 watchdog。先做语法检查、小范围 dry-run、BAM/VCF 可读性检查；失败只新建阶段目录，不清除任何旧文件。

## 7. 风险与停止条件

- 该重复区可能不能用短读长唯一定位；这是有价值的限制结论，不能以挑选 primary placement 掩盖。
- FISH 的片段观察不能直接换算为 WGS breakpoint；最终并列呈现而非合并计数。
- assembly 输入不足、BAM sequence 缺失比例高或 caller 源文件缺失时，量化缺失及对结论影响，继续可执行分支。
- 全窗口和旧 P41 candidate-gated 的事件 universe 若实质冲突，旧版不删除；新旧并列并解释差异来源。

## 8. 变更记录

- 2026-10-23 — 用户审阅完成并要求重新制定严谨计划、分析全部 18 样本；旧版 179 candidates/QNAME gate 及直接排除 secondary 的主线废止。旧结果/脚本复制归档；当前改为 full-window、molecule-aware、repeat-aware，新增 CNV、local assembly、DLD1 clone-specific 与大量图形交付。
- 2026-10-22 — 旧版建立；原文与实际脚本在顶部指定归档中保留，供逐项比较。
