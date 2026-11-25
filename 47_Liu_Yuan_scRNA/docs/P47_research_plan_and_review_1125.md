# P47 — Liu Yuan：293T 单细胞编辑效应分析

Plan Date（创建）: 2026-11-25

状态：数据下载和 MD5 校验已完成；GENEWIZ 公司分析结果已核对；本计划和唯一样本表已更新。原计划中“从本轮单细胞数据检出碱基编辑”的路线被测序化学否定。表达变化分析仍可能进行，但目前缺少可解释的样本设计；edited 与 control 标签又各自落在不同 GEX pool，暂不能做可归因的组间差异分析。正式分析待澄清目标和实验设计。

更新记录：

- 2026-11-25 — 按用户指示读取客户初始介绍并检查 GENEWIZ 项目 30-1371703644。下载日志及看门狗日志均确认本地 380/380 个文件、8/8 FASTQ MD5 通过、0 个失败。已将 8 个样本的 GENEWIZ probe barcode 与 Cell Ranger 细胞数写入 [唯一样本表](../scripts/0_sample_table.csv)。公司流程为 10x Fixed RNA Profiling（Flex）探针式 GEX，读取合成探针而非原始转录本碱基；更新分析决策表。未启动依赖样本设计的正式下游分析。

## 1. 客户目标、数据状态与来源

客户初始介绍指出 8 套 293T single-cell RNA-seq：4 套 mitochondrial gene-edited、2 套 base-editor controls、2 套 gRNA controls，目的是评估编辑效应。“mitochondrial gene-edited” 尚不能解释为 mtDNA 编辑，也可能是核编码的线粒体相关基因。客户所述类别来自客户介绍；下文将 GENEWIZ 样本标签暂按 edited/control 两类归组，具体 control 类型以及 N1/N2 的含义待确认。

### 1.1 下载与完整性

数据位于用户指定目录：

`/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/genewiz_30-1371703644_1124/30-1371703644/`

- 下载日志 `/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/logs/download_genewiz_20261124_0046.log` 于 2026-11-25 01:05:12 记录 `ALL_DONE`；8 个 FASTQ 的 `md5sum -c` 均为 `OK`。
- 看门狗日志 `/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/logs/download_watch.log` 于 2026-11-25 01:47:18 记录：本地 380 个文件、远端清单 380 个、MD5 `8/8`、失败 `0`、自动重启 `0`。
- 因此记录为“截至 2026-11-25 01:47:18，目录内 380 个文件”。
- 按远端清单，`00_fastq/` 有 16 个文件（8 个 FASTQ + 8 个 `.md5` sidecar，约 882.3 GiB）；`01_analysis/` 有 364 个文件（约 114.3 GiB）；全目录约 996.7 GiB。项目 `raw_data/` 中的 `genewiz_30-1371703644_1124` 符号链接指向上述原始目录。

### 1.2 公司分析流程和样本信息

两份 GENEWIZ `config.csv` 均记录 Cell Ranger 8.0.0、参考 `refdata-gex-GRCh38-2020-A`、探针面板 `Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv`，并以 Probe Barcode 拆分样本。交付目录包含 8 个样本的 raw/filtered feature-barcode matrices、molecule-info 文件和 QC 表；未发现 `.bam` 或 `.cram`。虽然 config 写有 `create-bam,true`，本地镜像并不含 BAM。

样本 ID、barcode、pool 与 `GEX: Cells` 来自 GENEWIZ config/QC。`group` 是按客户所述的 4 个 edited + 4 个 controls 暂定映射；控制类型未被公司 QC 表定义。样本表中的 `batch` 只记录公司 GEX pool 标签，不代表已确认的湿实验批次；target、biological replicate、clone 未提供，因此保留为 `unconfirmed`。

| 样本 ID | 暂定组别 | GENEWIZ GEX pool | Probe Barcode | GEX Cells | Median UMI/cell | Median genes/cell | Confidently mapped reads in cells |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| N1-ABE | edited | pool 1 | A-A10 | 15,245 | 47,228 | 8,625 | 97.61% |
| N2-ABE | edited | pool 1 | A-B10 | 19,027 | 52,216 | 8,867 | 97.99% |
| N1-CBE | edited | pool 1 | A-C10 | 20,706 | 48,541 | 8,878 | 97.70% |
| N2-CBE | edited | pool 1 | A-D10 | 15,726 | 58,424 | 9,175 | 97.75% |
| N1-Control | control | pool 2 | A-E10 | 20,780 | 61,946 | 9,109 | 97.83% |
| N2-Control | control | pool 2 | A-F10 | 20,053 | 61,535 | 9,199 | 97.80% |
| CBE-Control | control | pool 2 | A-G10 | 12,480 | 26,250 | 6,840 | 97.50% |
| ABE-Control | control | pool 2 | A-H10 | 12,027 | 42,071 | 8,209 | 97.87% |

公司报告合计 136,044 个 GEX cells，单样本为 12,027–20,780。公司过滤矩阵细胞数及 QC 指标是供应商结果，不是本项目独立重算值。FASTQ 按 pool 汇总；不能把同一个 pooled FASTQ 的大小或 read-pair 数填到多个样本行，故样本表中的 `r1_gib`、`r2_gib`、`raw_read_pairs` 留空。

## 2. 可行性和证据边界

### 2.1 当前 FRP 数据不能判定编辑碱基

这批数据是 10x Fixed RNA Profiling（Gene Expression Flex），不是测序读取原始 transcript-derived cDNA 的常规 3′/5′ 捕获。10x 官方说明：FRP 的 reads 来自与转录本杂交并连接的合成探针 DNA；若原转录本带有变异，该变异不会出现在测序数据中。因此，这批 FRP FASTQ、molecule-info 及可能存在的 probe-alignment BAM 都不能给出编辑位点 ref/alt allele 或编辑率。[10x：What is Fixed RNA Profiling?](https://www.10xgenomics.com/support/software/cell-ranger/7.2/getting-started/cr-flex-what-is-frp)；[10x：Fixed RNA Profiling with Cell Ranger multi](https://www.10xgenomics.com/support/software/cell-ranger/7.2/tutorials/running-pipelines/cr-flex-multi-frp)

现有数据最多支持探针目标基因的表达定量及探索性变化；不能用低表达、probe 捕获减少或表达矩阵变化替代 DNA/RNA 碱基级编辑证据。要回答编辑是否成功，需提供独立 DNA genotyping/amplicon 结果，或能够读取原始转录本序列的适配数据；本批 FRP 只适合表达分析。

### 2.2 表达组间比较存在 pool 混杂风险

**实测事实：** GENEWIZ 两份 config 分别列出 4 个 edited 标签和 4 个 control 标签样本，使用不同 GEX library/FASTQ pool。**推断：** 两个 GEX pool 可能对应不同 GEM/建库批次。如果确为不同 GEM/建库批次，组别与 pool 完全重合，组间差异无法与 pool 效应区分；统计校正不能从完全共线的设计中识别处理效应。需确认两 pool 的湿实验关系，及是否有跨 pool 的桥接/重复样本。

控制类别与独立样本结构也未确认。`N1-Control/N2-Control` 是否为 gRNA controls、`ABE-Control/CBE-Control` 是否为对应 editor-only controls，暂时仅从名称推测。N1/N2 可能是靶点、克隆或重复，当前不作假定。即便每类对照确为独立生物学样本，每类仍只有 2 个，低于本项目正式跨条件推断的 n≥3 独立样本门槛；最多谨慎呈现描述性结果。

### 2.3 不把 STARsolo 常规 3′/5′ 流程用于 Flex

现有 STARsolo 脚本及 10x v2/v3/v4 whitelist 验证针对常规 UMI 文库，不能直接用于 FRP 的 probe-based library。GENEWIZ 已按相应 probe set 用 Cell Ranger 8.0.0 完成样本拆分与表达定量。没有经验证的 FRP 二次计数方案前，不重跑 STARsolo，不把合成 probe reads 当作转录本序列，也不把公司基于 cluster 的 DE 结果当作样本级推断。

## 3. Decision table：本轮数据可支持的分析

| 研究问题 / data universe | 观察或证据单位 | 保留、去重与过滤规则 | 判定指标与 filter-impact audit | 当前允许的结论 |
| :--- | :---: | :---: | :---: | :---: |
| 编辑碱基是否出现；当前全部 FRP FASTQ 与 company outputs | 需要原始目标转录本碱基或 DNA molecule；当前读到的是 synthetic probe molecule | 不将 probe sequence 或 molecule-info 解释成 transcript sequence；不按表达阈值定义“已编辑细胞” | 标记为 assay-level `not assessable`；不生成 ref/alt、编辑率或“未编辑”计数 | 本批 FRP 不可用于碱基级编辑判定；结论不是“没有检测到编辑” |
| 每样本表达与文库 QC；8 份公司 raw/filtered matrices、Cell Ranger metrics、barcodes | sample × cell；跨条件单位为独立培养样本 | filtered matrix 表示公司 cell-called 集合；raw matrix 用于 cell-calling/ambient 背景核查；QC 按样本评估并记录每步 barcode 保留/排除原因 | 每样本 cells、UMI/基因分布、mapping/probe metrics；输出过滤前后 barcode 数及 excluded table | 可报告供应商 QC 与本项目的描述性表达 QC；不能推导编辑率 |
| edited vs control 的表达响应 | 独立生物学样本 × cell state × gene 的 raw-count pseudobulk；细胞数不能充当重复数 | control 类别保持分开；不合并两类 controls；不用 cell-level p 值伪增样本量 | 先审 biological replicate 与 GEX-pool 交叉；目前组标签与 pool 完全重合；报告 n、设计矩阵秩/共线性及过滤影响 | 当前不能归因组间差异；若设计确实完全混杂，只作逐样本描述，不做 inferential DE |
| 293T 中潜在状态异质性 | sample × cell；不是独立实验重复 | 逐样本评估 QC、doublet、ambient RNA；不把 UMAP 分离直接解释为细胞类型或编辑效应 | 报告每样本保留/排除数和原因；如用户确有该目的，再报告状态组成与敏感性 | 可做探索性细胞状态描述；不可泛化为组织细胞类型 |
| pathway / gene-set 变化 | 可识别对比下的 sample-level pseudobulk effect / ranked genes | 仅基于可解释的 contrast；不把公司 cluster-wise DE 当独立样本推断 | 报告背景集、multiple-testing、方向一致性与未通过的 contrasts | 只有设计可识别的表达对比才可进行探索性通路解释 |

后续每项过滤须给出前后 barcode/record 数量、原因及 machine-readable excluded audit。若某选择改变样本归组、观察单位或对比定义，先确认再运行。

## 4. 解除阻塞后的工作顺序

1. 补齐 8 个样本的真实组别、靶点、control 定义、独立培养/克隆/技术重复关系，以及两个 GEX pool 的实验关系；更新唯一 SSOT [0_sample_table.csv](../scripts/0_sample_table.csv)。
2. 确认本轮目标是表达响应，还是另有可读取原始 DNA/RNA 碱基的编辑验证数据。本批 FRP 不再安排碱基编辑检测分支。
3. 如做表达响应，从 Cell Ranger per-sample filtered/raw matrices 与 QC 起步，先做逐样本表达 QC、ambient 背景和 293T doublet 敏感性评估；不把常规 STARsolo 当作 FRP 重计数器。
4. 只有在独立样本数及 pool/batch 设计可识别时，建立 sample-level pseudobulk 对比；否则报告逐样本描述，不做伪重复显著性结论。
5. 按研究问题决定是否进行细胞状态与 pathway 分析；出图遵循本机 V21 palette 和科研作图规范，逐图核对样本点数、标签和可读性。
6. 完成英文研究报告及 `custom_research_report_YYYYMMDD/`，逐条审计主张、数据来源、过滤影响和结论边界；交付保持内部审阅、未发送。

## 5. 继续正式分析所需的信息

请补齐样本表中的信息；请勿仅依据 N1/N2 名称猜测：

1. **8 个样本的独立实验单位：** 每个样本是独立培养、独立克隆还是技术重复？N1/N2 分别代表靶点、克隆还是重复？
2. **4 个 edited 样本的构建信息：** 每个靶点的基因（mtDNA 或核编码）、参考 accession/坐标、ref 与预期 alt、编辑器和 gRNA；若目标实际是表达变化，也请说明。
3. **4 个 controls 的定义与对比：** N1-Control、N2-Control、ABE-Control、CBE-Control 分别包含什么构建体？哪些样本应互相比？
4. **两个 GEX pool 的实验定义：** 是否为不同 GEM、建库或处理批次？有无同一样本/参考样本跨 pool 作为桥接？
5. **结论目标：** 本轮是否只做 FRP 表达响应？若还要判断编辑是否成功，请提供已有 DNA amplicon/genotyping 数据的精确路径。FRP 本身不提供碱基编辑读数。

拿到设计信息后，继续运行可识别的表达分析；若设计仍完全混杂，则保留 QC 和逐样本描述，并在报告中明确哪些问题不可归因。

## 6. 当前完成状态

- 已完成：读取客户介绍；公司分析目录已在用户指定 raw-data 目录；日志确认本地 380/380、MD5 8/8 OK、0 FAILED；整理 8 个样本 barcode 与公司 QC；新增本计划。
- 尚未完成：样本设计澄清、可解释的表达差异分析、研究报告和交付包。正式表达分析未启动；这不是阴性分析结果。
- 不对外发送任何材料。后续如生成报告或交付包，默认仍待用户内部审阅。
