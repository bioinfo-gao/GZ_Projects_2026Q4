# P42 chr11 junction read-rescue：分析方法、结果与结论说明

- **文档性质**：内部中文技术说明；用于解释已经完成的 P42 `Analysis 1`，不是新的客户交付报告，也不替代逐候选结果表。
- **创建日期**：2026-10-22
- **分析范围**：S01–S10（明确包括后来追加的 S02、S03、S10）；不包括 Mutu I 或 DG75。
- **目标区间**：GRCh38 `chr11:114,604,000–114,625,500`。
- **本轮边界**：仅完成客户确认的 junction-level read rescue。未执行原信件中另行定义的 Analysis 2（S04/S05 与 S01/S02/S03 的克隆特异事件比较），也未执行 manual breakpoint curation、local assembly 或任何 wet-lab confirmation。

## 1. 一句话结论

在预先固定的候选筛选规则下，十个样本共有 179 个目标 junction 候选。将其原始支持 fragment 的序列重新提取后，以保留全部 alternative alignment 的方式独立重比对至 GRCh38，177/179 个候选能再次重建两个已记录 breakpoint 端：其中 4 个有至少 2 个独立 fragment 的重建支持，173 个只有 1 个 fragment 支持；2 个未重建。该结果说明大多数原始 read sequence 在相同参考基因组上具有可复现的 junction-level mapping 模式，但**不等同于确认 underlying structural variant（SV）真实存在**，尤其 173 个 single-fragment 结果不能单独承担 SV confirmation。

## 2. 本文回答的问题，以及不回答的问题

### 2.1 要回答的科学与技术问题

项目 41（P41）在 chr11 repeat-affected 区域产生了大量 junction 候选。短 read 在这类区域容易拥有多个相近的比对位置，因此仅凭一次 caller 或一次初始 split-read mapping，不能充分说明记录的两个 breakpoint 端能够由原始序列稳定地重现。

P42 Analysis 1 的问题被限定为：**对于每一个保留的 P41 junction，最初支持该 junction 的 fragment sequence 在独立重比对、且不隐藏 alternative alignment 的条件下，能否同时回到已记录的两个 breakpoint 端？**

这里的 `split read` 指同一个测序 fragment 的序列在参考基因组上以两个或多个相互分离的 alignment segment 表现出来；它可提示 junction，但本身不是独立的生物学重复。

### 2.2 本轮不回答的问题

以下问题超出本轮的预先授权范围，本文不会把它们写成已经得出的结论：

- 不判定任一候选是否为 confirmed SV，也不确定其 complete structure、allele fraction 或 clonality。
- 不比较 +dox 克隆与 DLD1 controls、也不比较不同 GM12878 clone 的 clone-specific event；这属于未执行的 Analysis 2。
- 不做 manual IGV review、base-pair-level breakpoint adjudication 或 local assembly。
- 不做 PCR 或其他 orthogonal assay；因此不产生 independent biological replication。
- 不根据样本名称推断客户未提供的实验元数据、细胞状态或表型。

## 3. 输入数据、坐标体系与可追溯性

### 3.1 输入来源

| 输入                                                       |                本轮用途                | 状态与来源                                                                             |
| :--------------------------------------------------------- | :-------------------------------------: | :------------------------------------------------------------------------------------- |
| 10 个 duplicate-marked BAM（Binary Alignment/Map）及其索引 |    提取原始支持 fragment 的 sequence    | `/home/gao/Dropbox/Vedita_Bam/`；十个 BAM 合计约 482 GiB；只读查询，未复制或修改 BAM |
| P41 target-evidence junction 表                            | 建立候选及记录初始 breakpoint / support | `/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs/analysis/region/evidence/target/`   |
| GATK GRCh38 analysis-set reference                         |             独立重比对参考             | P41 使用的同一 GRCh38 reference 与 `bwa-mem2` index                                  |
| P42 样本表                                                 | 本轮唯一的样本与预期候选数来源（SSOT） | `scripts/0_sample_table.csv`                                                         |

`samtools quickcheck` 已对 10 个纳入 BAM 逐一通过。P41 的已有交付目录保持冻结；P42 只在新建的 `analysis/20261022/` 与 supplementary delivery directory 中写入结果。

### 3.2 为什么本轮候选总数是 179，而不是旧邮件中的数字

本轮不直接沿用邮件中对八个样本的汇总数字，因为该数字没有附带可由程序重复执行的 eligibility rule。P42 在运行**前**将候选规则固定为：

1. P41 记录的 junction 必须标记为 `in_target == True`；
2. `n_split_nondup_mapq20 >= 1`，即至少有一个非 duplicate 的 split-read support 满足 P41 的两端 `MAPQ >=20` 口径。

再将这个规则应用于最终确认的 S01–S10 范围，得到 179 个候选。因样本范围已从原始八样本邮件扩展、且候选口径在运行前被明确化，两个汇总数不能视为同一分母下的直接差异。

### 3.3 候选数的样本分布

| 样本             |   客户标签   |    候选数    | Multi-fragment | Single-fragment | Not reconstructed |
| :--------------- | :----------: | :-----------: | :------------: | :-------------: | :---------------: |
| S01_DLD1_pop     |   DLD1_pop   |      20      |       1       |       19       |         0         |
| S02_DLD1_c1_ctrl | DLD1_c1_ctrl |      14      |       0       |       14       |         0         |
| S03_DLD1_c2_ctrl | DLD1_c2_ctrl |      13      |       0       |       13       |         0         |
| S04_DLD1_c3_dox  | DLD1_c3_dox |      17      |       0       |       17       |         0         |
| S05_DLD1_c4_dox  | DLD1_c4_dox |      12      |       1       |       10       |         1         |
| S06_GM12878_pop  | GM12878_pop |      22      |       2       |       19       |         1         |
| S07_GM12878_c5   |  GM12878_c5  |      24      |       0       |       24       |         0         |
| S08_GM12878_c6   |  GM12878_c6  |      13      |       0       |       13       |         0         |
| S09_GM12878_c7   |  GM12878_c7  |      26      |       0       |       26       |         0         |
| S10_GM12878_c8   |  GM12878_c8  |      18      |       0       |       18       |         0         |
| **总计**   |      —      | **179** |  **4**  |  **173**  |    **2**    |

所有样本的 `ambiguous_reconstruction` 和 `not_evaluable` 候选数均为 0。S02、S03、S10 在建立候选表阶段有显式 presence assertion，因此它们不会因早期八样本范围而被静默遗漏。

## 4. 分析设计与方法

### 4.1 总体设计

本流程不重新调用全基因组 SV，也不从新 reads 中寻找新的 junction。它是一个以 P41 已记录候选为对象的 read-rescue workflow：从原 BAM 回收产生该候选的原始 sequence，再把 sequence 放入一个独立的 alignment 环境中检验是否还能重建同一对 breakpoint 端。

工作流的逻辑如下：

```text
P41 target-evidence tables
  └─ 固定 eligibility rule → 179 个候选 junction
       └─ 按 P41 read name 从 BAM 扩展窗口回收原始 fragment sequence
            └─ bwa-mem2 -a -Y 对 GRCh38 重新比对
                 └─ 逐个候选检查两个已记录 endpoint（±20 bp）
                      └─ 汇总为 multi / single / ambiguous / not reconstructed / not evaluable
```

### 4.2 候选表的生成

脚本 `scripts/20261022/1_prepare_candidates.py` 逐样本读取 P41 的 `*.junctions.tsv`，仅保留第 3.2 节两项 eligibility criteria 同时满足的记录。生成的 canonical 表为：

- `analysis/20261022/candidates/candidate_junctions.tsv`：179 个候选、两端染色体/坐标/方向、P41 support counts 与 source read names；
- `analysis/20261022/candidates/candidate_support_fragments.tsv`：candidate-to-read-name 对照；
- `analysis/20261022/candidates/candidate_counts.tsv`：每个样本的候选数审计表。

候选类型的构成如下；类型来自 P41 原始表，并非 P42 新调用：

| P41 junction type  | 候选数 | 占 179 个候选的比例 |
| :----------------- | :----: | :-----------------: |
| BND（breakend）    |  155  |       86.59%       |
| INV（inversion）   |   20   |       11.17%       |
| DEL（deletion）    |   2   |        1.12%        |
| DUP（duplication） |   2   |        1.12%        |

### 4.3 原始 sequence 的回收

脚本 `scripts/20261022/2_extract_support_reads.py` 在扩展窗口 `chr11:114,554,000–114,675,500` 中查询每个 duplicate-marked BAM。该窗口比目标区间更宽，以避免仅在严格 target boundary 内查询而漏掉跨边界的支持 record。

具体处理为：

1. 只保留 read name 位于 P41 `source_reads` 清单的 alignment record；
2. 排除 secondary alignment；在查询窗口内若同一 read end 同时遇到 primary 与 supplementary record，则优先保留 primary record，以避免同一 read end 被重复写入 FASTQ；
3. 要求存在 query sequence； # 解释如下段
4. 将 sequence 与 quality 写为 per-sample FASTQ，并建立 `rescue_read_manifest.tsv` 将 FASTQ identifier 映射回 candidate、原 read name 和 mate。

#### “要求存在 query sequence”具体指什么

这里的 query sequence 是 BAM/SAM record 的 `SEQ` 字段：即该 read end 实际测得且随该 record 保存的碱基串；它不是参考基因组序列、不是该 read 的 CIGAR string，也不是“此 record 是否已经比对到某处”的判断。SAM/BAM 规范允许 `SEQ` 写为 `*`，其准确含义是**sequence 未被存储**；BAM 中对应的 `l_seq=0`。在 `pysam` 中，这表现为 `read.query_sequence` 为空或 `None`。P42 的下一步要把这个碱基串写成 FASTQ，再交给 `bwa-mem2` 从头搜索其可能的 reference placement；没有 `SEQ` 的 record 没有可供重新比对的输入，不能用参考碱基、CIGAR、同名的另一条 alignment 或人为补齐的 sequence 来替代。

常见情形与不应混淆的情形如下：

| 情形 | `SEQ` 是否可能为空 | 对本流程的含义 |
| :--- | :---: | :--- |
| 在 BAM 生成、转换、脱敏或节省存储空间时，有意丢弃 bases，仅保留 alignment metadata | 是 | 这是最直接的情形；必须回到保留原始 bases 的 FASTQ/BAM/CRAM 才可能重新比对 |
| 某些 alignment-only export、下游处理或不完整数据传输产生 record | 是 | record 仍可能有 QNAME、位置、CIGAR、FLAG 等信息，但不能成为 read-rescue 输入 |
| Unmapped record | 可以，但并非必然 | unmapped read 通常也可以保存原始 bases；`unmapped` 不等于 `SEQ` 缺失 |
| Secondary 或 supplementary alignment | 可以，但并非必然 | 这些是 alignment relationship 的 FLAG，不是 sequence-storage 标志；本流程因计数策略排除 secondary，但窗口内 supplementary record 可在没有对应 primary record 时提供 sequence |
| Low `MAPQ`、duplicate flag 或 vendor/QC flag | 否，二者没有必然关系 | 它们不是 `SEQ` 是否存在的判据；`MAPQ` 只在后续 unique reconstruction 判定时使用 |
| QUAL 缺失 | 否 | `QUAL='*'` 可以与非空 `SEQ` 共存；脚本会写入 quality placeholder，但不会补造 sequence |

另一个需要区分的情况是 **hard clipping（CIGAR 中的 `H`）**。hard-clipped record 的 `SEQ` 不是空的，但被 hard clip 的那部分 bases 不在该 record 的 `SEQ` 中；它与 `SEQ='*'` 的“完全未存储”不同。相反，soft clipping（`S`）的 clipped bases 仍在 `SEQ` 中。P42 的提取条件只要求 `SEQ` 非空，并未对每个 selected record 再作“是否含完整原始 read length”断言；因此流程优先选择查询窗口内的 primary record、仅在 primary 不在该窗口时采用可用 supplementary record，以尽量避免把 clipped representation 当作首选输入。该实现足以回答本轮已定义的 read-rescue 问题，但若某个候选需要 base-pair-level 结构判断，应进一步审阅对应 record 的 CIGAR、clipping 和原始 FASTQ，而不能仅以“`SEQ` 非空”作为完整 read 已恢复的证明。

这个筛选对应脚本中的条件 `if read.query_name not in wanted or read.is_secondary or not read.query_sequence: continue`。它依次排除三类 record：并非 P41 支持 read name 的 record、secondary alignment record、以及缺少实际 query bases 的 record。第三类的排除是一个**输入可用性条件**，不是对 read 质量、`MAPQ`、base quality 或 SV 真实性的过滤：

- `MAPQ` 在后续 endpoint reconstruction 阶段才用于判定 alignment 是否达到 `>=20`；本步骤不会因 `MAPQ` 低而删除 sequence。
- base quality 不会决定一个 `SEQ` 是否存在。若 BAM record 有 sequence 但没有可写出的 quality string，脚本会以长度相同的 `I` 作为 FASTQ quality placeholder；这只让 FASTQ 格式完整，不会补造或改变碱基序列。
- 对同一个 `(read name, mate)`，脚本在查询窗口内遇到 sequence 已存在的 primary 与 supplementary record 时优先保存 primary record。这一选择的作用是每个 read end 只写一次 FASTQ，避免同一 read end 由多个 alignment record 重复计数。脚本不会拼接多个 alignment record 的 sequence，也不会从 primary / supplementary record 之一推断另一条 record 缺失的碱基；实际写入的是被选中 record 的 `SEQ` 字段。若 primary record 本身不落在查询窗口内，窗口内的 supplementary record 仍可作为该 read end 的可用 sequence 来源。

若一个候选的所有预期 source fragment 都因缺少 `SEQ` 而无法写出 FASTQ，该候选最终会是 `not_evaluable`，而不是 `not_reconstructed`：前者表示没有可以测试的原始序列，后者表示序列已被测试但不能同时重建两个端点。实际运行中，179 个候选均至少回收到一个 source fragment，185 个按原始 read name 计的 source fragment 均进入了 re-extraction 计数，因此没有 candidate-level `not_evaluable` 结果。脚本没有单独输出被跳过的单条 BAM record 清单，故这个结果只能说明没有候选因缺少可用 sequence 而失去可评估性；不应额外表述为“每一条落在查询窗口内的 BAM record 都含 sequence”。

179 个候选的 `source_fragments` 总数为 185，且 185/185 均在 BAM 中重新提取到，因此没有 `not_evaluable` 候选。因为 paired-end fragment 可能贡献两个 read end，FASTQ manifest 有 298 行；候选层面的 fragment 数始终按去重后的原始 read name 计数，避免把同一 fragment 的两个 end 误作两个独立分子。

### 4.4 独立重比对

每个样本的 FASTQ 以 `bwa-mem2 2.2.1 mem -a -Y -t 8` 比对到 P41 所用的 GATK GRCh38 analysis-set reference。

| 参数     |                       实际作用                       | 本项目中的必要性                                                              |
| :------- | :--------------------------------------------------: | :---------------------------------------------------------------------------- |
| `mem`  |          执行 short-read sequence alignment          | 与 P41 的短读长参考比对框架保持可比                                           |
| `-a`   | 输出所有找到的 alignment，包括 alternative alignment | repeat-affected region 中不能因只保留一个 best alignment 而掩盖多重比对       |
| `-Y`   |   supplementary alignment 使用 soft clipping 表示   | 保留 split/supplementary alignment 的 query-sequence 表达，便于 endpoint 检查 |
| `-t 8` |                  使用 8 个计算线程                  | 仅影响运行资源与速度，不改变判定阈值                                          |

各样本 SAM 在保留第一个 header 后合并为 `analysis/20261022/remap/all_samples.sam`，供统一的逐候选分类脚本读取；分类仍由 manifest 中的 sample/candidate identifier 限定，样本之间不会混合计数。

### 4.5 Endpoint matching 与分类规则

脚本 `scripts/20261022/3_classify_rescue.py` 对每一个 FASTQ read end 的全部重比对记录检查 P41 记录的两个 breakpoint endpoint。

对于候选的每一端，必须同时满足：

1. alignment reference chromosome 与记录的 `chr1` 或 `chr2` 相同；
2. 使用与 P41 `dir1` / `dir2` 一致的 alignment endpoint：`+` 端比较 `reference_end`，`-` 端比较 `reference_start + 1`；
3. 该 endpoint 与 P41 记录位置的绝对差不超过 20 bp；
4. 用于左右端的两个 alignment record 必须不同，避免把同一条普通 alignment 同时冒充两个断点端。

若有多个能够命中两端的 alignment pair，代码选择 alignment score（`AS` tag）之和最高的 pair。若这对 alignment 的两个 `MAPQ` 都至少为 20，则该 read name 形成 `unique_reconstruction`；若可以命中两个 endpoint、但任一 `MAPQ <20`，则为 `ambiguous_reconstruction`。这里的 `MAPQ` 是 mapping quality，而不是 coverage 或 fragment count。

最终 candidate-level status 的优先级和定义如下：

| Candidate status              |                                    规则                                    | 对证据强度的含义                                               |
| :---------------------------- | :-------------------------------------------------------------------------: | :------------------------------------------------------------- |
| `multi_fragment_supported`  |            至少 2 个不同原始 read name 为 unique reconstruction            | 同一 candidate 有多于一个独立测序 fragment 的 read-level 支持  |
| `single_fragment_supported` |            恰好 1 个不同原始 read name 为 unique reconstruction            | 同一 candidate 只由一个 fragment 重建；不能视作 confirmed SV   |
| `ambiguous_reconstruction`  | 无 unique reconstruction，但至少一个 fragment 可重建且不满足 MAPQ criterion | 端点模式可见，但 mapping ambiguity 未被阈值排除                |
| `not_reconstructed`         |       已提取 sequence，但没有任何可同时匹配两个端点的 alignment pair       | 本轮 read rescue 不支持该候选；不自动等价于 biological absence |
| `not_evaluable`             |                  没有从 BAM 回收到该候选的 source fragment                  | 输入 sequence 缺失，无法进行 read-rescue 判定                  |

一个 candidate 若同时含 unique 与 ambiguous fragment，candidate-level 结果以 unique reconstruction 为准；ambiguous fragment 数仍被完整保留在结果表。这一优先级解释了本轮“0 个 ambiguous candidate”与“5 个 ambiguous reconstructed fragment”可以同时成立：5 个 ambiguity fragment 都属于另有至少一个 unique fragment 的 candidate，没有候选仅停留在 ambiguous 层级。

## 5. 结果

### 5.1 总体 read-rescue 结果

| Candidate-level 结果                     |    候选数    |       比例       | 解释                                    |
| :--------------------------------------- | :-----------: | :--------------: | :-------------------------------------- |
| `multi_fragment_supported`             |       4       |      2.23%      | 至少 2 个不同 fragment 独立重建两端     |
| `single_fragment_supported`            |      173      |      96.65%      | 仅 1 个 fragment 独立重建两端           |
| `ambiguous_reconstruction`             |       0       |      0.00%      | 无仅具 ambiguous mapping 的候选         |
| `not_reconstructed`                    |       2       |      1.12%      | 提取到原 sequence，但未成对重建两个端点 |
| `not_evaluable`                        |       0       |      0.00%      | 所有候选均有可评估的 source fragment    |
| **至少一个 unique reconstruction** | **177** | **98.88%** | read-rescue 在本预设阈值下重建两个端点  |

在 fragment 层面，185 个重新提取 fragment 中有 182 个产生 unique reconstruction；另有 5 个 fragment 出现了不满足 `MAPQ >=20` 的 endpoint reconstruction。后者不是 5 个额外 candidate，因为这些 fragment 分别属于已有 unique reconstruction 的 candidate；所有 fragment-level 数目应以逐候选表的字段定义为准，不能将不同字段跨层级相加。

### 5.2 Multi-fragment-supported 候选

下表列出证据更强的 4 个候选。它们是“较强的 read-level 支持”，不是自动获得 SV confirmation 的事件。

| Candidate            | Type | 已记录 breakpoint                                     | Source / re-extracted / unique fragments | 解释                           |
| :------------------- | :--: | :---------------------------------------------------- | :--------------------------------------: | :----------------------------- |
| S01_DLD1_pop_J001    | BND | `chr2:69,105,893 (-)` ↔ `chr11:114,609,953 (+)`  |                2 / 2 / 2                | 两个原始 fragment 均按阈值重建 |
| S05_DLD1_c4_dox_J001 | BND | `chr11:114,607,192 (+)` ↔ `chr21:20,648,241 (-)` |                3 / 3 / 3                | 三个原始 fragment 均按阈值重建 |
| S06_GM12878_pop_J001 | BND | `chr2:42,949,036 (+)` ↔ `chr11:114,606,321 (-)`  |                2 / 2 / 2                | 两个原始 fragment 均按阈值重建 |
| S06_GM12878_pop_J002 | BND | `chr3:104,481,114 (-)` ↔ `chr11:114,617,659 (-)` |                2 / 2 / 2                | 两个原始 fragment 均按阈值重建 |

### 5.3 未重建候选

| Candidate            | Type | 已记录 breakpoint                                      | Source / re-extracted / unique fragments | 本轮判定              |
| :------------------- | :--: | :----------------------------------------------------- | :--------------------------------------: | :-------------------- |
| S05_DLD1_c4_dox_J010 | DEL | `chr11:114,614,044 (+)` ↔ `chr11:114,622,342 (-)` |                1 / 1 / 0                | `not_reconstructed` |
| S06_GM12878_pop_J018 | BND | `chr2:123,920,435 (+)` ↔ `chr11:114,611,172 (+)`  |                1 / 1 / 0                | `not_reconstructed` |

两个候选的 source sequence 均成功从 BAM 回收，因此它们不是数据缺失；在本轮 reference、全 alternative alignment 输出、端点方向与 ±20 bp tolerance 下，均没有形成一个可同时命中两端的 alignment pair。这个结果只应表述为“本轮未重建”，不应升级为“候选一定不存在”或“P41 一定错误”。

### 5.4 出现 ambiguous fragment 的候选

虽然没有 candidate 被最终归类为 `ambiguous_reconstruction`，以下 5 个 single-fragment-supported candidate 同时保留了 1 个 ambiguous reconstructed fragment。该信息提示其 alternative alignment landscape 仍需在任何高风险生物学决策前审阅：

- `S05_DLD1_c4_dox_J002`
- `S06_GM12878_pop_J012`
- `S07_GM12878_c5_J010`
- `S08_GM12878_c6_J001`
- `S09_GM12878_c7_J005`

它们各自仍有至少一个 `MAPQ >=20` 的 unique reconstruction，故最终 status 保持 `single_fragment_supported`。逐 read / 逐 fragment 的替代比对关系未被简化为“完全唯一”；`-a` 的目的正是让此类竞争 alignment 不被输出策略隐藏。

## 6. 质量控制与验证

本轮执行的验证包括：

| 验证项               | 实际检查                                                                    | 结果                                          |
| :------------------- | :-------------------------------------------------------------------------- | :-------------------------------------------- |
| BAM 完整性           | 对十个 in-scope BAM 执行 `samtools quickcheck`                            | 10/10 通过                                    |
| Scope 完整性         | 脚本断言 S02、S03、S10 必须在候选表中出现                                   | 通过                                          |
| 候选数               | 脚本断言最终候选行数为 179                                                  | 通过                                          |
| Source fragment 回收 | 每个候选至少回收到一个 source fragment                                      | 179/179 通过；185/185 source fragments 已回收 |
| 分类全覆盖           | 逐候选 status 汇总必须覆盖全部 179 行                                       | 4 + 173 + 0 + 2 + 0 = 179                     |
| 脚本完整性           | Python scripts 执行 `py_compile`；shell scripts 执行 `bash -n`          | 通过                                          |
| 正式交付一致性       | analysis result 与 supplementary delivery 中复制的 candidate table 逐行比较 | 179 行完全一致                                |
| 文件格式             | 正式输出检查为 LF line ending，无 CRLF                                      | 通过                                          |

结果文件的 SHA-256 可用于日后核对是否发生非预期变化：

- `analysis/20261022/results/rescue_candidate_results.tsv`：`a04be987f742a333a14ffa600f74aae19ac1a364c68923dc41c1b1d7c2638c60`
- `analysis/20261022/results/rescue_sample_summary.tsv`：`7f2a224ca7974fd085a5d7266aacc57c942d42b1e062131ab386ea9bcb081b4`

## 7. 如何正确解读这些结果

### 7.1 “重建成功”实际证明了什么

`multi_fragment_supported` 或 `single_fragment_supported` 证明的是：在所定义的 reference、alignment program、alternative-alignment output、方向和 ±20 bp endpoint tolerance 下，原 P41 source fragment 的 sequence 还能产生同时符合两个记录端点的 alignment。它加强了“这条 sequence 可产生该 mapping pattern”的证据。

它**没有**独立证明该 breakpoint 存在于足够多的细胞或等位基因中，也没有确定两个端点之间的全长结构。因为 P42 重用了 P41 的原始 sequence，这是一种 computational re-evaluation，而不是一个新的 biological replicate。

### 7.2 Multi-fragment 与 single-fragment 的实际区别

四个 multi-fragment-supported 候选至少由两个不同 read name 重建，因而比 single-fragment observation 更不容易完全由单个随机 read 或某一个 alignment record 驱动。它们应是后续人工审阅或 targeted confirmation 的优先候选。

173 个 single-fragment-supported 候选构成本轮绝大多数。它们的重比对能够重复到相同 endpoint，但 evidence unit 仍只有一个原始 fragment。对于 repeat-affected region，这种证据足以保留为 read-level observation，却不适合直接转换为“confirmed rearrangement”“clonal event”或样本间生物学差异的结论。

### 7.3 为什么 0 个 ambiguous candidate 不能被误读成“没有 mapping ambiguity”

candidate-level `ambiguous_reconstruction` 仅指“没有任何 unique reconstruction、但至少一个 fragment 可重建且 MAPQ 不达阈值”。本轮 0 个候选落在这一专门分类，说明每个可重建候选至少有一个满足 `MAPQ >=20` 的 endpoint pair。

这并不代表所有 read 都只有一个可能位置。五个候选仍出现了 ambiguous reconstructed fragment，且整个设计特意使用 `-a` 来输出 alternative alignment。因此正确结论是“最终 candidate label 未被 ambiguity 单独主导”，不是“该 repeat-affected region 没有多重比对风险”。

### 7.4 两个未重建候选应如何处理

S05_DLD1_c4_dox_J010 与 S06_GM12878_pop_J018 不应列入“read-rescue-supported”候选。由于两个 fragment 均已实际提取，这个结论不是由 input absence 造成；但短 read sequence、repeat context、alignment scoring 与 endpoint tolerance 均可能影响是否能重建固定端点。没有进行 base-pair curation 或 assembly 的情况下，不应以这两个结果反推事件的绝对 biological absence。

### 7.5 与 P41 的关系

P42 是对 P41 chr11 目标区 junction 证据的一层补充与细化：它把已有候选的 source sequence 放入明确定义的 alternative-alignment-aware rescue test。P41 的冻结交付没有被修改。本轮增加的是每个候选的可重建性证据和证据强弱分层，而不是以没有执行的工作替代或否定此前记录。

## 8. 局限性与后续使用边界

1. **同源数据局限**：P42 重用 P41 的 source reads；重新比对可以检验 mapping reproducibility，不能替代独立测序、独立建库或独立样本。
2. **short-read 局限**：单个短 fragment 无法在此流程中恢复完整 rearrangement architecture，也无法可靠判断 distant partner 与目标区的完整连接关系。
3. **repeat context 局限**：尽管 `-a` 保留 alternative alignment，`MAPQ >=20` 是预定义的 operational threshold，并不把所有重复相关 ambiguity 消除。
4. **endpoint tolerance 局限**：±20 bp 是本次固定的判定窗口；它适合评估与 P41 记录端点的一致性，不等同于 base-pair-resolved breakpoint 定义。
5. **定量局限**：本流程不估计 variant allele fraction、细胞比例、copy number 或群体/克隆频率。
6. **比较设计局限**：本轮不进行 DLD1 +dox vs control、或 GM12878 clone 间的统计比较，因此不能从本结果推出 treatment effect 或 clone specificity。

若某个候选将被用于生物学决策，建议以 4 个 multi-fragment-supported 候选为优先级起点，再进行 candidate-specific manual alignment review 与跨越候选 junction 的 targeted orthogonal confirmation。上述工作在本轮未执行，产生新结果时应写入新的、带日期的 supplementary analysis directory，而不是覆盖现有交付。

## 9. 文件索引与复现入口

| 文件（相对于项目根目录）                                   | 角色                                                                 |
| :--------------------------------------------------------- | :------------------------------------------------------------------- |
| `scripts/0_sample_table.csv`                             | 样本范围、BAM 路径和每样本预期候选数的 SSOT                          |
| `scripts/20261022/0_run_all.sh`                          | 完整正式 workflow 的 orchestrator；在 `tmux` 内运行                |
| `scripts/20261022/1_prepare_candidates.py`               | 从 P41 frozen target-evidence 建立 179 行候选表                      |
| `scripts/20261022/2_extract_support_reads.py`            | 从十个 BAM 的扩展窗口回收 source sequence                            |
| `scripts/20261022/3_classify_rescue.py`                  | 基于全部重比对记录执行 endpoint matching 和 candidate classification |
| `analysis/20261022/candidates/candidate_junctions.tsv`   | 候选及原始 P41 support 信息                                          |
| `analysis/20261022/reads/rescue_read_manifest.tsv`       | FASTQ identifier 与 candidate/read/mate 的映射                       |
| `analysis/20261022/results/rescue_candidate_results.tsv` | **机器可读的逐候选 canonical result**                          |
| `analysis/20261022/results/rescue_sample_summary.tsv`    | 每样本汇总                                                           |
| `custom_research_report_supplementary_20261022/`         | 已生成、待用户过目的英文客户交付包；尚未发送                         |

本说明中的汇总数字均来自 `analysis/20261022/results/rescue_candidate_results.tsv` 与 `rescue_sample_summary.tsv`。如需追踪任意一句结论，应先在这两张 canonical table 中定位 candidate ID 或 sample ID，再回到 P41 frozen evidence table 和对应 BAM 验证原始记录。
