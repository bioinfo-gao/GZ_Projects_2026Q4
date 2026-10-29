# P47 — Liu Yuan：293T 单细胞编辑证据分析

Plan Date（创建）: 2026-10-29

更新记录：
- 2026-10-29 — 阅读 Gmail、调研 GitHub 原仓库、建立待确认的证据模型与配色方案。仅调研，未启动正式分析。

## 1. 已知事实与输入状态

来源：[Gmail：Quote of Single-Cell RNA-seq Data Analysis](https://mail.google.com/mail/u/?authuser=zorin.gao%40gmail.com#all/1a18816b9908228d)，2026-10-29；其中引用客户 2026-10-28 原信。
客户说明：8 套 293T single-cell RNA sequencing（scRNA-seq）数据，其中 mitochondrial gene-edited 4 套、base editor control 2 套、guide RNA（gRNA）control 2 套。核心问题：能否从这些数据检测编辑效应。

邮件中的分析承诺是工作目标，不是已证明的数据可行性。尚不能把 mitochondrial gene-edited 等同于已确认编辑 mitochondrial DNA（mtDNA）：也可能指核编码的线粒体相关基因，必须先拿到确切位点。

本轮仅检查既有项目目录与配色来源，尚未定位/清点客户原始数据。每样本文件大小、read pairs、细胞数及总数据量均未测量，不能填 0。待定位送样单时，以同目录作为 raw data 路径，记录“截至日期时间目录内文件数”，再逐样本汇总。样本唯一来源为 `../scripts/0_sample_table.csv`；当前仅建表头，不虚构客户样本 ID 或重复结构。

## 2. GitHub 现有方案与初步取舍

| 方案与原仓库 | 适用角色 | 本项目初步选择 |
| :--- | :---: | :---: |
| [maegatk](https://github.com/caleblareau/maegatk) | 面向 unique molecular identifier（UMI）scRNA 数据的线粒体变异计数与分子 consensus，特别支持 MAESTER | mtDNA 靶点成立且 BAM 标签适配时优先评估；普通 scRNA 覆盖能力需实测 |
| [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite) | 从 BAM 按 barcode/UMI 做位点 allele pileup，支持指定候选及 discovery | 已知靶位点的通用计数候选，也适合核编码靶点；与主方法对照不构成独立实验验证 |
| [mgatk](https://github.com/caleblareau/mgatk) | 主要面向线粒体单细胞基因分型，尤其 ATAC 数据 | 作为参考，不直接把 ATAC 默认过滤和灵敏度搬到 RNA |
| [MQuad](https://github.com/single-cell-genetics/MQuad) | 从线粒体变异识别克隆信息 | 非主线；已知编辑位点不需要先通过克隆信息筛选才计数 |
| [Scanpy](https://github.com/scverse/scanpy) | 表达质控、降维、聚类 | 表达主线候选，保留 raw counts 与全部基因 |
| [muscat](https://github.com/HelenaLC/muscat)、[pertpy](https://github.com/scverse/pertpy) | 多样本 differential state / perturbation 分析框架 | 借鉴 sample-level 聚合与扰动比较，不机械套用 Perturb-seq guide assignment |

以上为 2026-10-29 在线文档审阅；未安装、未做本机兼容性或性能验证。没有单一仓库自动解决实验设计、位点可检测性与结论边界。

## 3. 拟议 decision table：正式运行前确认

| 研究问题 / data universe | 观察与证据单位 | 保留与去重规则及理由 | 指标与 audit | 允许的结论 |
| :--- | :---: | :---: | :---: | :---: |
| 每个预期靶点是否可检测；从每样本全部输入 alignment 出发，覆盖所有预定靶点，不以旧变异表封闭输入 | sample × barcode × locus × molecule；read/alignment 为来源证据 | 保留 primary/secondary/supplementary、FLAG、坐标、CIGAR、MAPQ、base quality、CB/UB、read-end；barcode/UMI/位置及冲突处理规则需按 assay 适配，不能仅按 UMI 全局合并 | 全 barcode 与表达 QC 后 barcode 分别统计 coverage；各过滤前后 record/read-name/molecule/candidate 数；低质量、歧义、无标签、分子冲突逐项列 audit | 可评估 / 覆盖不足；无覆盖是 NA，不是未编辑 |
| 可评估位点是否出现预期碱基转换 | 独立分子计数；样本才是跨培养重复单位 | PCR duplicates 不多计分子；原始 alignment context 留存。核基因组与 mtDNA 联合检查定位，关注 nuclear mitochondrial sequences（NUMTs）、链/位置偏倚；不能把 MAPQ 或 primary 当真来源证明 | ref/alt/other 分子、RNA allele fraction、有效覆盖、对照背景、置信区间；预期转换方向按参考链归一；阈值需覆盖审计后预注册 | RNA 层面与预期编辑一致的证据；不直接等同 DNA 编辑率或所有细胞的编辑比例 |
| 是否有转录状态变化；全部合格表达细胞，不限有靶位点覆盖者 | 单细胞用于状态描述；sample × state raw-count pseudobulk 用于样本比较 | 逐样本自适应质量控制（QC）、doublet/ambient RNA 评估；不得仅因线粒体比例高删除可能有生物学信号的细胞；保留过滤原因。表达分支的 ambient 校正不替代 allele 分支的分子证据 | 每样本保留/排除细胞与原因；未整合基线；batch 校正敏感性；细胞周期与应激状态；不默认回归掉线粒体比例 | 293T 状态分布与通路关联；不硬套免疫细胞图谱、不把 UMAP 分离视为编辑证明 |
| edited 对两类 control 的表达差异 | 独立培养样本为重复，不以细胞数当 n | 两类 control 先分开；先确认 4 edited 是同靶点重复还是不同编辑。按本地 SOP，每组独立 n≥3 才正式跨条件推断；pseudobulk 每格≥10 细胞并列出被排除格 | 设计矩阵、batch/condition 共线检查；当前两类 control 各仅2套，按描述性效应量与一致性报告，不出正式显著性结论 | 可描述变化；不能用合并 control 绕过设计门槛 |
| 编辑负荷是否与表达相关 | 有足够靶点覆盖的细胞，保留样本归属 | 单独分析覆盖选择偏倚、UMI/测序深度、细胞周期；不把未覆盖细胞标成 unedited | covered 与 uncovered 的表达/QC 分布；分样本趋势与敏感性 | 探索性关联，不是因果或独立生物学验证 |
| 非靶向变化 | 只有经确认范围内、RNA 实际覆盖的碱基位点 | 如扩展 discovery，须另行定义候选宇宙、错误模型、多重检验和完整过滤审计 | 可评估区域与未覆盖区域分开 | 表达扰动不等于 off-target mutation；scRNA 不能排除全基因组脱靶 |

SEQ 非空与完整原始 read 分别审计：记录 hard/soft clipping、序列长度、read-end；需要时回查 FASTQ。所有重比对仅验证 mapping/algorithmic reproducibility。

## 4. 执行顺序（待输入与科学口径明确）

1. 收齐样本设计与靶点；清点 FASTQ / Cell Ranger BAM、索引、raw/filtered matrices、barcodes、reference、web_summary。只有表达矩阵则先报告无法做碱基证据分析。
2. 优先做 8 套数据的靶位点覆盖与 barcode/UMI 标签审计，决定能否做单细胞、仅样本汇总或不可评估；不先跑完整表达工作流才发现位点无覆盖。
3. 并行设计 allele 与表达两条分支；各自保留过滤 audit，通过 sample + barcode 连接。
4. 先交可检测性结果，再给 RNA 编辑证据、状态变化和限制。无数据前不承诺耗时/灵敏度。

## 5. 两套配色的采用方式

唯一色值来源：`/home/gao/Utilities/color_and_format/figure/palettes.csv`。不复制色表；Python/R 使用现成 `palette_priority.select_palette()`。

- Categorical_7：三个实验臂暂固定 named mapping：edited、base_editor_control、gRNA_control，按 `select_palette(3, scheme="Categorical_7")` 取色；得到真实设计后一次固化，不逐图重新取色。状态类别的独立语义另建映射并做跨 panel 撞色检查。
- Paired_4：用于真实配对设计，或以显式 `scheme="Paired_4"` 的 independent 模式取最多4类颜色。当前不能因恰有8套数据就强行画成4对。设计未确认前保留这套方案为有条件采用，不伪造配对。
- 8 个样本先用实验臂颜色 + sample 标签/形状，不硬挤进7分类色表。小点/细线采用 `role="stroke"`，填充/大点用 fill；文字黑色，箱线 median 黑色。未覆盖单元格用灰色/NA；allele fraction 用独立连续色标，不用 categorical 彩虹。
- 固定变量跨图同色；语义不同 panel 颜色不相交；样本级图核对8点是否齐全。实际出图时再做白底可读性与逐 panel 视觉检查，本轮未生成生物学图。

## 6. 正式分析前需要确认的关键输入

1. 8 个真实样本 ID、独立培养/克隆/技术重复关系；4 个 edited 的靶点、预期转换、分组与 batch；两类对照分别包含什么构建体及应比较关系。
2. 靶点是 mtDNA 还是核编码基因，参考 accession/坐标/参考与预期碱基；编辑器、guide 或其他定位机制、处理时间。客户原文 gRNA 的具体作用尚未定义，不据此推断编辑器类型。
3. 可用数据是否包含 BAM/FASTQ 与 barcode/UMI 信息，是否有可关联的已有 DNA 编辑检测结果。没有 DNA 结果也可评估 RNA，但结论层级须保持。

当前状态：调研完成、设计待输入；未运行正式程序、未安装软件、未发送邮件。
