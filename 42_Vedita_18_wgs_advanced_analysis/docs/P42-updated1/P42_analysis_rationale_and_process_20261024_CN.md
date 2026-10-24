# P42 分析 rationale 与过程说明（清晰版）

**建立日期：** 2026-10-24  
**用途：** 内部研究、复核与交接。它重新说明 P42 的问题、证据模型、执行过程和结论边界；不替代任何历史记录，也不修改既有结果。  
**原始计划归档：** `docs/P42/P42_research_plan_and_review_1022.md`  
**当前客户问题英文报告：** `custom_research_report_supplementary_20261024_v8/P42_Client_Request_Completion_1024.md`

## 1. 一句话问题定义

P42 的核心不是泛化地“找出所有 SV”，而是在 18 个 short-read WGS 样本中，以 full-window、repeat-aware 的方式复核 P41 报告第 7.3 节的 chr11 repeat-region 线索，判断是否存在可被严格定义为 clone-/doxycycline-specific 的 junction，并明确哪些证据只能支持“候选假设”、不能支持“已确认 breakpoint”。

客户随后特别关心：P41 的 EV0623 在 P42 中是什么、C2 的 junction 是否更明确、Manta/DELLY 是否一致、unmapped-mate rescue 是否增加支持、DLD1 control 与 +doxycycline 样本有何差异，以及可否交付 BAM/caller data。P42 的后续工作因此是对 P41 §7.3 的深化，不是另一项独立的 discovery study。

## 2. 必须先区分的四个层次

| 层次 | 观察单位 | 不能混同为 | 在 P42 中的用途 |
| :--- | :---: | :---: | :--- |
| Alignment record | 一条 primary、secondary 或 supplementary record | 独立 molecule | 保留 CIGAR、FLAG、MAPQ、soft clip、SA 与 mate context |
| Read end | `(sample, QNAME, read end)` | 独立 fragment | 描述 split alignment、query overlap/gap 与 read-level support |
| Fragment/template | `(sample, QNAME)` | Independent biological replicate | 作为 candidate support 的主要计数单位；raw 与 standard-eligible 分开 |
| Biological event | 多条可定位、相互一致、满足敏感度要求的 fragments，并可有独立 assay 支持 | 单一 best-hit 或单个 caller record | 只有通过全部预设 gates 才能称为 high-confidence event |

这一区分是 P42 最重要的设计选择。Secondary/supplementary records 不能作为新增独立 molecules 叠加计数，但它们也不能因“避免双计数”而被删除，因为它们是 repeat、split/chimeric structure 与 placement ambiguity 的直接证据。

## 3. 数据 universe 与固定比较对象

### 3.1 分析 universe

- 18 个 GRCh38-aligned WGS BAM；mean depth 约 8.2x–19.3x。
- chr11 repeat-region full-window 的原始 alignment context，而非只从旧 caller table 或旧 candidate list 开始。
- 保留 primary、secondary、supplementary、MAPQ 0、multi-mapping、soft-clipped、unmapped mate；任何过滤均须有 before/after audit。
- Local assembly、contig multi-placement、Manta、DELLY、GRIDSS/SvABA 的 outputs 作为互补 evidence streams，而不是互为 biological replicates。

### 3.2 固定 candidate universe

完整 evidence inventory 产生 282 个 candidate observations，其中 94 个具有可唯一解析的双端 geometry。针对 DLD1 S01-S05，四个 case-only candidates（C1-C4）进入固定 re-genotyping/comparison。固定集合的目的，是让所有 DLD1 样本接受同一坐标、同一支持规则和同一 depth context 的比较；它不允许把未进入该集合的全窗口 evidence 反过来写成“不存在”。

## 4. 研究问题与预注册式判定表

| 研究问题 | 证据与保留原则 | 数值/结构判定 | 能说什么 | 不能说什么 |
| :--- | :---: | :---: | :--- |
| 是否有 repeat-aware mapping evidence？ | 从 full-window raw records 开始；保留全部 alignment context 和 filter-impact audit | Evidence clusters、candidate geometry、assembly placement | 某位置有 mapping/junction-like evidence | 单个 primary placement 就是真实 locus |
| P41 EV0623 是否对应 P42 C2？ | 对照双方 chr11/partner geometry、sample、support 与 representation | P41: chr11:114,607,192–chr21:20,648,241；P42: chr11:114,607,191;chr21:20,648,242 | 两端各 1 bp 的表示差异下，是同一 geometry 的 carry-forward | 这是两个独立事件 |
| C2 的 read-end representation 是否一致？ | 所有可用 sequence-bearing records，按 QNAME/read end 保留 | 5 个 read ends、3 QNAME、MAPQ 60、`AT` 2-bp overlap、0-bp gap | 与 2-bp microhomology 一致的 split-read representation | 唯一、经生物学验证的 breakpoint |
| Manta/DELLY 是否 concordant？ | 72 个 source inputs 不丢记录 | 两 BND ends 在任意方向均须落在 candidate 两端 ±50 bp | 0/20 exact Manta、0/20 exact DELLY；C4 chr11-only anchors 非 concordance | “没有 exact record”就是生物学不存在 |
| 是否 clone-/doxycycline-specific？ | S01-S05 采用同一 fixed-candidate rule；raw/standard counts 分开 | ≥3 standard-eligible QNAME、≥60x local context、sequence-level junction | 只有三项全过才为 high-confidence | 在 14.7x–23.1x control depth 的 zero counts 下宣布高灵敏度阴性 |
| Rescue 是否改变结论？ | 5 个 rescued unmapped mates 与每个 candidate QNAME set 相交 | C1 1→1、C2 3→3、C3 2→2、C4 1→1；intersection 均 0 | 本次 rescue 对 C1-C4 的增量为 0 | Rescue 在整个 locus 没有价值 |

## 5. 实际执行过程

### 阶段 A：全窗口 evidence inventory 与 coverage context

首先建立 18 样本的全窗口 record/molecule/mate-rescue inventory，并同时报告不同 FLAG、MAPQ、sequence availability 与 filter impact。随后以分层 coverage/CNV tracks 观察 chr11 window 的 depth context。此阶段的产物回答“输入中看到了什么”，而不是提前给 event 定性。

### 阶段 B：候选聚类、assembly 与 placement ambiguity

原始 junction-like evidence 被聚类；local assembly 只对有信息量或 caller-anchored 的 bins 执行，但未选 bins 始终保留在 audit 中。Contigs 以全 placement 模式比对，记录 alternative placements；因此 contig 不能仅凭一个 best hit 被提升为唯一插入位点。完成后共有 59 assemblies、90 contigs、32 alternative placements。

### 阶段 C：Caller-compatible evidence 与 fixed-candidate comparison

Manta/DELLY source records、GRIDSS/SvABA context、assembly/read evidence 被统一对照，但每一类证据的独立性边界被保留。对 C1-C4，S01-S05 使用同一两端 geometry 和同一 counting rule 做固定 re-genotyping，避免 case/control 因不同筛选而不可比较。

### 阶段 D：P41 §7.3 专项深化

1. 将 P41 EV0623 与 P42 C2 crosswalk，确认不是新事件；
2. 在 C2 的 5 个 read ends 上逐一计算 query overlap/gap；
3. 以 two-end 50-bp rule 审计 Manta/DELLY；
4. 以 candidate-QNAME intersection 衡量 rescue 改变量；
5. 将 support、standard eligibility、depth 与 sequence-level evidence 一起判定 DLD1 clone/doxycycline gate；
6. 仅进行 C2 in-silico PCR design 与 specificity screen，不把设计当成 PCR result。

### 阶段 E：交付整合

当前正将主 P42、assay follow-up 和客户问题深化结果合并到单一新目录 `custom_research_report_20261024/`。该阶段不重算科学结果：它复制 18 BAM/BAI、标准化并脱敏 caller VCF headers、将 DELLY BCF 转为 indexed VCF、汇集 14 张正式图和 machine-readable tables，并进行 checksum、format、header-path 与交付边界审计。历史交付目录保留不改。

## 6. 当前结果应如何解读

| 对象 | 观察到的事实 | 当前可用结论 | 关键限制 |
| :--- | :---: | :---: | :--- |
| 全队列 | 282 candidate observations、94 个双端几何可解析、4 个 DLD1 case-only candidates | 没有 candidate 被升格为 high-confidence SV | Repeat placement ambiguity、短读长 coverage context |
| C1 | S05 raw/standard support 1/1 | Low-support candidate | 无独立 junction/clone-specific proof |
| C2 | S05 raw 3，standard 0；5 read ends 均 `AT` 2-bp overlap/0-bp gap | P41 EV0623 的优先 carry-forward hypothesis | 无 standard-eligible support、深度不足、无独立 confirmation |
| C3 | S04/S05 raw 1/1，standard 0/1 | 两个 +doxycycline samples 的 low-support candidate | 未过 high-confidence gate |
| C4 | S05 raw/standard 1/1；5 条 DELLY chr11-only anchors | Low-support candidate，anchor context 被保留 | chr11-chr3 BND 未得 exact concordance |
| Manta/DELLY | 20 comparisons 的 exact two-end match 都为 0 | 无 caller exact concordance | 同一 WGS 的 caller comparison 不是独立 validation |
| Rescue | 5 recovered mates，所有 C1-C4 intersections 为 0 | 对这四个 candidates 无增加 | 不外推到其他 loci |

## 7. 已建立、未建立与下一步

**已建立：** P41 EV0623→P42 C2 的正确 crosswalk；C2 的 read-end microhomology-like representation；caller exact-concordance 的缺失；rescue 对 C1-C4 的零增量；四个 fixed candidates 都未过预设 high-confidence gate；C2 的 in-silico assay design。

**未建立：** 唯一真实 genomic breakpoint、clone-/doxycycline-specific biological event、经验证 allele fraction、与 FISH/cytogenetic subpopulation 的因果联系、project-wide SNV/indel callset 的存在。

**合理的后续证据：** 在 independent DLD1 genomic DNA 上测定 C2 junction-spanning PCR assay，并以相应 controls 和实测 analytical sensitivity 解释；若希望达到 3-fragment/5% target，先确认 local depth 是否达到 60x context；若要连接 FISH，须取得 raw images 与 per-cell scoring table。这里没有也不应以重新运行同一短读长数据取代独立验证。

## 8. 关键文件导航

| 用途 | 当前文件 |
| :--- | :--- |
| 客户问题英文报告 | `custom_research_report_supplementary_20261024_v8/P42_Client_Request_Completion_1024.md` |
| 中文阅读版 | `docs/P42-updated1/P42_Client_Request_Completion_1024_CN.md` |
| C2 read-end 证据 | `custom_research_report_supplementary_20261024_v8/results/c2_read_end_junction_resolution.tsv` |
| Caller comparison | `custom_research_report_supplementary_20261024_v8/results/candidate_manta_delly_concordance.tsv` |
| DLD1 gate comparison | `custom_research_report_supplementary_20261024_v8/results/clone_specific_comparison_with_criteria.tsv` |
| Rescue accounting | `custom_research_report_supplementary_20261024_v8/results/read_rescue_candidate_impact.tsv` |
| C2 assay design | `custom_research_report_supplementary_20261024_v5/results/assay_design_selection.tsv` |
| 单一客户包（构建中） | `custom_research_report_20261024/` |

本说明只重组 rationale 与过程；数值主张均可回溯到上述 source tables 或正式英文报告。若有新的原始 assay/FISH/SNV-indel input，应作为新阶段输入，不覆盖本轮 evidence ledger 或结论边界。
