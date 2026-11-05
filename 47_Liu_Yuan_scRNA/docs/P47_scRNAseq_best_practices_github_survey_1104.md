# 单细胞转录组分析当前最佳实践：GitHub 生态调研

调研日期：2026-11-04
定位：内部方法学文档，供 P47 与后续单细胞转录组（single-cell RNA sequencing，scRNA-seq）项目选型。
状态：只做了调研和玩具数据验证。未读取客户数据，未运行正式分析，未安装或升级软件。

本文不重复已有 SOP。静默错误清单、环境路径和断言库见 [`/scrna` skill](/home/gao/.claude/skills/scrna/SKILL.md)，P47 的证据模型见 [P47_research_plan_and_review_1029.md](P47_research_plan_and_review_1029.md) §3。本文只写选型依据、2026 年的新证据和这些证据带来的改动。

## 0. 结论

**常规多样本项目的推荐路线如下；P47 不能照搬，差异见 §6。**

| 环节 | 默认选择 | 替代或补充 |
| :--- | :---: | :---: |
| 定量 | Cell Ranger；需要开源流程时用 nf-core/scrnaseq 4.2.0 | STARsolo、alevin-fry、kallisto-bustools |
| 环境 RNA（ambient RNA） | 污染明确：CellBender 或 DecontX full mode；污染不明或只有 filtered 矩阵：SoupX reduced mode | 校正结果另存一层，不覆盖原始 counts |
| 双细胞（doublet） | 逐样本 scDblFinder | 多样本混样时加 hashing 或基因型拆分 |
| 质量控制（quality control，QC） | 逐样本中位数绝对偏差（median absolute deviation，MAD）标记，先审计再过滤 | miQC 联合建模线粒体比例与文库复杂度 |
| 归一化与特征选择 | shifted log；2,000 个高变基因（highly variable genes，HVG）起步，只打标记 | analytic Pearson residuals、deviance |
| 整合 | 先看未整合基线。简单批次用 Harmony 或 Seurat，复杂任务用 scVI 或 scANVI | scib-metrics 量化，并单独查过度校正 |
| 跨条件差异表达（differential expression，DE） | 样本 × 细胞类型 raw counts 求和后用 pseudobulk + DESeq2 或 edgeR | 重复测量或随机效应用 dreamlet；muscat 可批量封装 |
| 细胞组成 | sccomp 或 scCODA | speckle（propeller）；无标签时用 Milo |
| 富集 | decoupler、fgsea | 过表达分析（over-representation analysis，ORA）按 UP/DOWN 分开做 |

**这轮调研带来五个需要改变做法的发现：**

1. **scDblFinder 在同质细胞群里的能力接近只看文库大小。** 本机模拟中，细胞类型异质时 recall 为 0.78–0.93，同质时只有 0.41–0.44（§5）。这对 293T 这类单一细胞系有直接影响。
2. **本机 scDblFinder 1.24.10 与 xgboost 3.4.2.1 的组合落在上游修复之前。** 上游 1.27.9 记录了 xgboost ≥3.4 问题的修复，并取消了训练失败时的静默回退。本机阳性对照能正常识别，但正式项目仍应先升级或在日志里确认未回退（§5）。
3. **环境 RNA 校正没有全能工具。** 2026 年的三项 benchmark 一致认为不存在通用最优；scAR 和 CellClear 被指会重构 count 矩阵；低污染数据上 CellBender 和 DecontX full mode 更容易过度校正（§3、§4）。
4. **Seurat 与 Scanpy 的结果差异已经正式发表。** 两者默认设置的差异，约相当于只测不到 5% 的 reads 或只分析不到 20% 的细胞所带来的变异（检索摘要，引用前核对原文）。一个项目应固定一条主线，并记录版本（§4）。
5. **“深度学习扰动预测赢不过线性基线”已不是定论。** 2026 年 Nature Biotechnology 论文认为换用校准后的指标结论会变。我们仍不把预测结果写进结论，但理由应改为“有争议且未独立验证”（§4）。

## 1. 调研范围与证据等级

**调研方法：** 通过 GitHub API 查询 44 个候选仓库名，其中 3 个发生迁移、2 个 scrapper 候选名不存在；另做主题检索，并阅读关键 README、NEWS、changelog 和 sc-best-practices 的 notebook 原文。元数据原表见 [github_repo_metadata_1104.tsv](P47_scRNAseq_best_practices_evidence_1104/github_repo_metadata_1104.tsv)。

**Star 数只反映知名度，不代表方法正确。** 本文的选型依据是 benchmark、维护状态和本机实测。

| 标记 | 含义 | 使用要求 |
| :--- | :---: | :---: |
| 【实测】 | 在本机运行得到，日志在证据目录 | 可直接复核 |
| 【原文】 | 直接读取 GitHub 上的 README、NEWS、changelog 或 notebook | 可直接引用 |
| 【页面】 | 通过网页抓取工具读取论文页面后摘要，未逐字通读 | 关键数字引用前核对 |
| 【检索】 | 只见到搜索结果摘要，未打开全文 | **对外引用前必须读原文** |
| 【判断】 | 我的方法学意见 | 可以反对 |

**调研缺口：**

- 一条关于 293T 与线粒体 DNA（mitochondrial DNA，mtDNA）编辑的检索被模型安全策略拦截，我没有改写措辞重试。因此本文**不包含**293T 核型、亚克隆结构或该类编辑器的文献结论。
- Nature Communications（Nie 等）和 Cell Genomics 的全文未能打开（需登录或返回 403），相关内容均标为【检索】。
- 未检索到 2025–2026 年独立、同行评议的 doublet 工具头对头 benchmark，也未检索到细胞系线粒体 QC 的专门 benchmark。

## 2. 开跑前的证据模型

全局规则 §2A 要求每项分析先写清证据模型。下面是常规项目的模板，P47 以自己的计划文档为准。

| 研究问题 | 证据单位 | 保留与去重 | 判定指标 | 允许的结论 |
| :--- | :---: | :---: | :---: | :---: |
| 有哪些细胞状态或类型 | 细胞 | 保留 flagged 细胞及原因，不直接丢弃 | 簇稳定性、marker、注释置信度 | 描述性状态；簇的存在本身不是独立证据 |
| 跨条件表达是否变化 | 样本 × 细胞类型 pseudobulk | 求和 raw counts；<10 细胞的格子列出后剔除 | 效应量、错误发现率（false discovery rate，FDR）、样本间一致性 | 仅在独立样本足够时才给统计结论 |
| 组成是否变化 | 样本内比例 | 保留所有注释类别，含 Unassigned | 成分模型的置信区间 | 比例变化；不直接推断增殖或死亡 |
| 批次校正是否有效 | 细胞嵌入 | 与未整合结果并排 | 批次指标与生物保留指标分别报告 | 校正后可比较；不等于差异已被证明 |
| 位点或等位基因证据 | sample × barcode × locus × molecule | 分子去重用于计数，原始 alignment 保留作来源证据 | 有效覆盖、ref/alt 分子数、对照背景 | 仅限可评估位点；无覆盖是 NA，不是阴性 |

## 3. 逐环节最佳实践

### 3.1 上游

- **Cell Ranger：** 从 v7 起默认计入内含子 reads，`--include-introns` 默认为 true；`--create-bam` 是必填参数【检索】。需要做位点或等位基因分析时必须设为 `true`。
- **参考基因组：** 同一项目的全部样本用同一个参考包。10x 当前人类参考 GRCh38-2024-A 使用 GENCODE v44【检索】。
- **流程化：** nf-core/scrnaseq 4.2.0 支持 simpleaf、STARsolo、kallisto + bustools 和 Cell Ranger【原文】。nf-core/scdownstream 的设计值得参考，但 GitHub Releases 上未见正式版本，暂不用于交付【原文】。
- **RNA velocity：** 需要 spliced 和 unspliced 分层。velocyto.py 的最新 release 还是 2019 年，优先考虑 alevin-fry 或 kb-python【判断】。

### 3.2 环境 RNA

| 场景 | 选择 | 依据 |
| :--- | :---: | :---: |
| 明确污染，且有 raw 矩阵和 GPU | CellBender | 综合排名第一；需要 GPU，运行时间长，内存高【页面】 |
| 明确污染，无可用 GPU | DecontX full mode | 综合排名第二；去除较少，但利于生物学解释【页面】 |
| 污染不明，或只有 filtered 矩阵 | SoupX reduced mode | 低污染或阴性对照上比前两者更不易过度校正【页面】 |
| 不建议作默认 | scAR、CellClear | 预印本报告 CellClear 替换 >93% counts，scAR 产生未校正数据中不存在的细胞类型【页面】 |

**做法：**【判断】把校正结果另存为 layer，原始 counts 不动，并报告校正前后的差异。校正不一定改善数据（Cargnelli 2026 明确提到）。不同 benchmark 对 scAR 的评价有分歧，见 §4。

本机注意：GTX 1060 3 GB 不足以支撑 CellBender 的常规使用（沿用 skill 的判定）；CellBender 0.4.0 仅更新 Python 兼容性，输出与 0.3.x 一致【原文】。

### 3.3 Doublet

- **默认：** scDblFinder，逐样本运行。书中引用的 Xi & Li 补充分析中其准确度最高【原文】。
- **需要知道的局限：**
  - README 明确说明同类型细胞形成的 homotypic doublet 在转录组上几乎无法区分，通常无法用该包识别【原文】。
  - F1000Research 的评审者指出其参数在同一批 benchmark 数据集上调过【检索】。
  - Cell Genomics 2024 用合成 DNA 条形码做真值，发现数量接近正确，但 AUPRC（area under the precision-recall curve）偏低【检索】。
- **做法：**【判断】把 doublet 分数当连续风险，不当事实标签；文库大小上尾单独做敏感性分析；混样项目优先用 hashing 或基因型信息。

### 3.4 QC

- **sc-best-practices 的具体写法【原文】：** `log1p_total_counts`、`log1p_n_genes_by_counts`、`pct_counts_in_top_20_genes` 用双侧 5 MAD；线粒体比例经 logit 变换后取上侧 3 MAD，且同时超过 8%。书中说明 8% 并非通用值。
- **线粒体比例要结合其他指标看。** 书中强调三项 QC 指标必须联合判断【原文】。miQC 以 75% 为默认风格的丢弃概率，联合建模线粒体比例和文库复杂度【检索】；它并不总比固定阈值宽松。
- **做法：**【判断】先写 `qc_flag_*` 布尔列，看过 flagged 细胞在哪些簇、哪些样本后再决定过滤。这样过滤原因可追溯，也符合 §2A 的 filter-impact audit。

### 3.5 归一化、特征选择与降维

- **归一化：** shifted log + 主成分分析（principal component analysis，PCA）仍是稳妥默认。Ahlmann-Eltze & Huber 2023 比较 22 种变换，没有其他方法能持续胜出【检索】。analytic Pearson residuals 适合选择生物学可变基因和稀有细胞类型【原文】。
- **特征数：** Zappia 2025 用 2,000 个特征，理由是它在各数据集和指标上都稳定得分较高，同时建议按用途调整；query 任务需要更多特征【页面】。
- **batch-aware HVG：** 同一研究未发现一致优势，不必默认使用【页面】。
- **警惕“混得很好”：** 只含噪音的特征集也能得到很好的批次混合【页面】。因此批次指标高不代表整合成功。
- **scale：** Luecken 2022 发现 scaling 会让方法偏向去批次、牺牲生物保留【检索】。

### 3.6 整合

- **选型：** 书的结论是不存在全场景最优。Harmony 和 Seurat 适合简单批次校正，scVI、scGen、scANVI、Scanorama 适合复杂任务【原文】。
- **差异很大的系统**（跨物种、类器官与组织、细胞与细胞核）：考虑 sysVI，并通过 cycle-consistency 权重调节强度；该方法也存在过度整合的风险【检索】。
- **评估：** 用 scib-metrics 分别报告批次去除和生物保留。标准 scIB 指标对细胞类型内部变异不敏感，Genome Biology 2025 的 scIB-E 专门扩展了这一点【检索】。
- **边界：**【判断】整合后的嵌入只用于聚类和可视化。DE 始终回到原始 counts。

### 3.7 聚类与注释

- **聚类：** Leiden，`n_iterations=2`；结果写入带分辨率的 key；不解读 UMAP（uniform manifold approximation and projection）上的簇间距离【原文】。Scanpy 1.12 已弃用 `sc.tl.louvain`【原文】。
- **过度聚类：** 先聚类再在同一批数据上做 DE 属于 double dipping，同质群体也能被切开并“找到” marker【检索】。需要声称“这是真实亚群”时，加一层校准：
  - count splitting（Neufeld 2024），但过离散越强，效果越接近 double dipping【检索】；
  - recall 或 ClusterDE，二者的 benchmark 均由作者本人完成【检索】。
- **注释：** 手动 marker、分类器（CellTypist、SingleR）和参考映射三条线都走。书中建议 LLM（large language model，大语言模型）和基础模型最多作为额外的第一轮筛选【原文】。

### 3.8 跨条件 DE

- **默认：** pseudobulk。书中明确写到样本层面方法优于细胞层面【页面】；Squair 2021、Junttila 2022 也支持这一点【检索】。Scanpy 1.12 文档已加入警告，指出 `rank_genes_groups` 不适用于大多数场景【原文】。
- **具体做法【页面】：**
  - 格子阈值的常用指引是 ≥10 细胞且 ≥1,000 总 counts，书中说明阈值与数据集相关；
  - 每个细胞类型分别过滤低表达基因；
  - 先对 pseudobulk 做 PCA，把与主要成分相关的协变量放进设计矩阵。
- **随机效应：** dreamlet 已于 2026 年发表于 Nature Communications，用精度加权线性混合模型处理 pseudobulk【检索】；开发者自评，尚未见独立复核。本机未安装。
- **Python 还是 R：** pertpy 中的 edgeR 仍通过 rpy2 调用 R；PyDESeq2 作者自述结果“相似但不完全相同”【检索】。最终结论沿用 R 原版，与 skill 一致。

### 3.9 组成、富集、通讯、轨迹与调控网络

| 分析 | 选择 | 分寸 |
| :--- | :---: | :---: |
| 组成 | sccomp、scCODA、propeller、Milo | 书中提醒重复数少时置信区间很宽，且条件与混杂共线时 k 近邻（k-nearest neighbor，kNN）方法会丢信息【原文】 |
| 富集 | decoupler、fgsea | 连续统计量优先用 t 值而非 log2FC【原文】 |
| 细胞通讯 | LIANA+、CellChat | 结果是表达相容性假说；sc-best-practices 0.2.0 新增了专章【原文】 |
| 轨迹 | CellRank、Palantir 等 | scVelo、CellRank 均维护中；velocity 依赖上游分层 |
| 调控网络 | pySCENIC | 最新 release 为 2022 年，最近一次 push 为 2025-06；依赖钉死，存在维护风险【实测】 |
| 扰动筛选 | pertpy（Mixscape、E-distance） | 书中的流程面向有 guide assignment 的 pooled screen【原文】 |

## 4. 证据之间的分歧

| 问题 | 一方 | 另一方 | 本文取舍 |
| :--- | :---: | :---: | :---: |
| pseudobulk 是否总是最优 | Squair 2021、Junttila 2022：控制假阳性更好【检索】 | Nguyen 2023：批次效应大时 pseudobulk 表现差；GLIMES 2025：过于保守（作者自有方法）【检索】 | 仍作默认。阴性结果写成未检出，不写成无变化 |
| scAR 是否可用 | Nie 2026、Cargnelli 2026：去除能力强【检索/页面】 | Schroeder 2026 预印本：产生假细胞类型；Cargnelli 也报告其损失内源信号【页面】 | 不作默认 |
| DecontX 的评价 | Nie 2026：综合分高，细胞间变异保留好【检索】 | 另有来源指出其对高污染基因校正不足【检索】 | 高污染基因单独复查 |
| Seurat 与 Scanpy 能否互换 | 相同输入矩阵时 log 归一化结果等价【检索】 | Rich 2026：HVG、图构建、p 值并列处理不同，结果差异明显【检索】 | 不混用、不互相当验证，版本写入 Methods |
| 深度学习扰动预测 | Ahlmann-Eltze 2025：未持续胜过线性基线【检索】 | Nature Biotechnology 2026：校准指标下可以胜出【检索】 | 有争议，不进结论句 |

## 5. 本机差距与实测

### 5.1 版本差距

| 组件 | 本机 | 上游最新 | 影响 |
| :--- | :---: | :---: | :---: |
| scanpy | 1.11.5（Python 3.11.16） | 1.12.5 | 1.12 已移除 Python ≤3.11，**不能原地升级，需新建环境** 【原文/实测】 |
| scvi-tools | 1.3.3 | 1.5.2 | 项目内固定版本即可 |
| Seurat | 5.5.1 | 5.6.0 | 5.6.0 修复 `PercentageFeatureSet`，本机实测未受影响（见下） |
| scDblFinder | 1.24.10（xgboost 3.4.2.1） | devel 1.27.9 | **有风险，见下** |
| anndata / pertpy / decoupler | 0.12.19 / 1.0.3 / 2.2.0 | 0.13.4 / 1.5.0 / 2.3.1 | 跨大版本 API 可能变化，升级前先核对 |

**`sc_tut_r` 和 `sc_tut_main` 中未安装【实测】：** SoupX、celda（DecontX）、sccomp、dreamlet（以上在 R 环境中查询）；scib-metrics（Python 环境中查询）。skill 把 SoupX 和 DecontX 列为环境 RNA 的默认工具，但 R 全栈环境 `sc_tut_r` 里没有。`DE_R45` 本轮未查。

**未找到 Cell Ranger：** 在 `PATH`、两个环境的 `bin`，以及 `/home/gao/Utilities`、`/opt`、`/mnt/ex_8T_SSD` 下 3 层以内都没找到。这不代表其他位置没有。本机 STAR 为 2.7.10a。

### 5.2 Seurat 线粒体比例

Seurat 5.6.0 的 changelog 记录了 `PercentageFeatureSet` 读取 layer 的修复。【实测】本机 5.5.1 在默认 assay、`pattern="^MT-"` 下，未拆分和按样本拆分的对象都与手算真值完全一致（最大差异为 0）。这只覆盖这种常规调用，未覆盖 `assay=`、`layer=` 等其他参数组合。

### 5.3 scDblFinder 阳性对照

上游 NEWS 的 1.27.9 记录包括修复 xgboost ≥3.4 问题、取消训练失败时的静默回退、修复 doublet 来源标签错误【原文】。本机是旧版加新 xgboost，所以我做了已知真值的模拟。

【实测】模拟 3,000 细胞 × 2,000 基因，真 doublet 由两个细胞的 counts 相加得到，默认参数，每种条件 2 个种子：

| 条件 | 真 doublet 比例 | recall | precision | AUPRC | 只用文库大小的 AUPRC |
| :--- | :---: | :---: | :---: | :---: | :---: |
| 3 个细胞类型，只保留跨类型 doublet | 6.8–7.1% | 0.78–0.93 | 0.92–0.96 | 0.95–0.98 | 0.76 |
| 单一同质群体 | 11.1% | 0.41–0.44 | 0.91–0.92 | 0.85–0.86 | 0.83–0.84 |

**解读：**

- **本机组合能工作。** 跨类型条件下明显优于文库大小基线，没有训练完全失败的迹象。但模拟无法排除真实数据上出现回退，因为旧版不报告回退。
- **同质群体中，它几乎只用到了文库大小信息。** AUPRC 只比文库大小基线高 0.01–0.04。对单一细胞系，应把漏检写进 Limitations。
- **这是偏乐观的玩具模拟。** 没有真实的细胞周期、深度结构或 ambient，且 doublet 由严格相加生成；Cell Genomics 2024 指出这种生成方式有利于该方法【检索】。数字只说明方向，不能当作真实数据的检出率。

**`doublet 来源` 列在本机版本中不可信：** README 说明 1.27.8 之前该标签大多是错的【原文】。不要使用这一列。

### 5.4 pseudobulk 文件桥

【实测】在玩具数据上跑通 Python 导出到 R DESeq2 的完整链路，并断言 pseudobulk counts 总和等于细胞级 counts 总和。这只验证代码链路，不验证任何生物学结论。脚本：[check1_pseudobulk_export.py](P47_scRNAseq_best_practices_evidence_1104/check1_pseudobulk_export.py)、[check2_pseudobulk_deseq2.R](P47_scRNAseq_best_practices_evidence_1104/check2_pseudobulk_deseq2.R)。

```python
pb = sc.get.aggregate(adata, by=["sample", "state"], func="sum", layer="counts")
# 守恒断言：聚合不应丢失或凭空产生 count
assert pb.layers["sum"].sum() == adata.layers["counts"].sum()
# scanpy 1.11 没有 n_obs_aggregated（1.12 才加入），需要自己按 sample × state 计数后合并到 obs
```

## 6. 对 P47 的适配

P47 是 293T 单一细胞系，8 套数据（edited 4、base editor control 2、gRNA control 2），靶点和分组尚待客户确认。证据模型仍以[计划文档](P47_research_plan_and_review_1029.md) §3 为准，下面只写常规路线需要改动的地方。除标明来源外，均为【判断】。

| 常规做法 | P47 的处理 | 原因 |
| :--- | :---: | :---: |
| 细胞类型注释 | 不做图谱映射，只描述细胞周期、应激和线粒体相关状态 | 细胞系映射到组织图谱没有意义 |
| 默认做批次整合 | 主线用未整合基线；Harmony 只作可视化敏感性分析 | 单一细胞系的主要差异可能就是处理效应，整合有抹除风险 |
| scDblFinder 过滤 | 仍运行，但结合文库大小上尾报告，并在 Limitations 写明同质群体的漏检 | §5.3 |
| 环境 RNA 校正进入 DE | 默认不校正，只估计污染水平并做敏感性分析 | 单一培养物中 ambient 谱与细胞谱接近，校正可辨识性低；低污染数据上工具有过度校正风险【页面】 |
| 线粒体比例过滤 | 只标记不直接过滤；过滤前后各做一次对比 | 该比例可能正是处理效应本身；书中的 8% 也声明不通用【原文】 |
| 回归线粒体比例和细胞周期 | 默认不回归，只记录并做敏感性分析 | 可能回归掉待检信号 |
| 同一数据上聚类再报 marker | 声称存在亚群前先做校准（count splitting 或 ClusterDE） | 同质细胞系容易被切出假亚群；ClusterDE 在 5 个同质细胞系上做过测试（作者自评）【检索】 |
| 跨条件统计 | 对照两类各 2 套，低于 n≥3 门槛，按计划只做描述性效应量 | 还需确认 8 套是独立培养还是技术重复 |
| Mixscape 区分响应细胞 | 不作主线 | 它依赖 guide assignment；按表达选响应细胞再检验表达有循环风险 |

**如果靶点确认是 mtDNA，等位基因分支有几条硬约束：**

- **未富集的 3′ 数据覆盖很低。** MAESTER 原始研究中，普通 scRNA-seq 的单细胞 mtDNA 平均覆盖约 0.2–0.7×，富集后为 52–217×【检索】。
- **覆盖随位点变化很大。** 一项 2026 年的成纤维细胞研究中，MT-ND5 的平均覆盖为 3.11–5.61，MT-ND4 为 34.90；tRNA 和 rRNA 不被捕获【页面】。因此可评估性必须逐位点实测。
- **样本层面比单细胞层面可靠。** 同一研究中，scRNA pseudobulk 与 DNA 异质性相差 0–6 个百分点，但低覆盖会高估 0% 和 100% 的极端值【页面】。建议先做样本级等位基因比例，再评估能否下沉到单细胞。
- **不要照搬该研究的阈值。** 它要求 ≥100 条线粒体 reads、位点 ≥2 reads，但它的测序深度为每细胞 16–23 万 reads，且针对已知变异【页面】。P47 的阈值要在覆盖审计后预先写定。
- **maegatk 的适配性要先验证。** 它默认丢弃 <3 reads 的 UMI（unique molecular identifier）以降低 PCR 错误【检索】，在低深度数据上可能损失大量分子；仓库最近一次 push 是 2024-01【实测】。cellsnp-lite 仍在维护（2026-09）【实测】。两者结果对照只验证算法可重复性，不是独立生物学验证。
- **可加一个 ambient 背景对照。** 对 raw 矩阵中的空液滴同样做靶位点 pileup，估计 ambient 分子带来的等位基因背景，而不是只比较对照样本。

## 7. 未采纳的项目

| 项目 | 未采纳原因 |
| :--- | :---: |
| omicverse、SCP | 一体化封装，默认值和版本隐藏在封装层；鉴于版本差异已被证明影响结果，交付流程不走这层【判断】 |
| scGPT 等基础模型 | 书中定位为最多第一轮筛选；2026 年基础模型 benchmark 仍以 scVI、Harmony 作基线【原文/检索】 |
| mLLMCelltype | 注释工具，P47 不做细胞类型注释；常规项目也只作第一轮 |
| medical-research-skills、SciAgent-Skills 等 agent skill 合集 | Star 高但属于提示词合集，不是方法学证据 |
| rapids-singlecell | 本机 GPU 为 GTX 1060 3 GB，不具备条件 |

## 8. 待决定事项

以下都没有执行。

| 事项 | 建议 | 回退方式 |
| :--- | :---: | :---: |
| 为 P47 准备 R 环境：安装 SoupX、celda，升级 scDblFinder | 正式分析前做。用新环境，不动 `sc_tut_r`，后者是已验证的教程环境 | 删除新环境 |
| 新建 Python 3.12 环境以使用 scanpy 1.12 | 暂不需要，P47 期间固定 1.11.5 并写入 Methods | 无需回退 |
| 更新 `/scrna` skill | 建议改四处：scDblFinder “benchmark 最强”改为带局限说明；补 homotypic 与同质群体判据；环境 RNA 选型按场景分；扰动预测改为“有争议”。**未修改，等你确认** | git 回退 |
| 更新 P47 计划 §2 的 maegatk 行 | 补上维护状态和低深度风险。正式运行前修订 | git 回退 |

**P47 启动仍缺客户输入：** 8 个样本的真实 ID 与独立性、靶点坐标（mtDNA 或核基因）、对照定义、是否有 BAM 或 FASTQ。这些见计划 §6。本轮未定位原始数据，样本表仍只有表头。

## 9. 主要来源

**GitHub 与官方文档**

- [theislab/single-cell-best-practices](https://github.com/theislab/single-cell-best-practices)（v0.2.0，2026-10-29，未经正式同行评审）及[在线版](https://www.sc-best-practices.org/)
- [scverse/scanpy 1.12 release notes](https://github.com/scverse/scanpy/tree/main/docs/release-notes)、[satijalab/seurat NEWS](https://github.com/satijalab/seurat/blob/main/NEWS.md)
- [plger/scDblFinder README](https://github.com/plger/scDblFinder) 与 [NEWS](https://github.com/plger/scDblFinder/blob/devel/inst/NEWS)
- [nf-core/scrnaseq](https://github.com/nf-core/scrnaseq)、[nf-core/scdownstream](https://github.com/nf-core/scdownstream)
- [broadinstitute/CellBender](https://github.com/broadinstitute/CellBender)、[caleblareau/maegatk](https://github.com/caleblareau/maegatk)、[single-cell-genetics/cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite)
- [10x Cell Ranger release notes](https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-release-notes)、[参考包说明](https://www.10xgenomics.com/support/software/cell-ranger/10.0/release-notes/cr-reference-release-notes)

**论文与预印本**

- 总论：[Heumos 2023, Nat Rev Genet](https://www.nature.com/articles/s41576-023-00586-w)
- 版本差异：[Rich 2026, Cell Systems](https://www.cell.com/cell-systems/abstract/S2405-4712(26)00042-6)（[预印本](https://www.biorxiv.org/content/10.1101/2024.04.04.588111v1)）
- 环境 RNA：[Cargnelli 2026, Front Bioinform](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2026.1844838/full)、[Nie 2026, Nat Commun](https://www.nature.com/articles/s41467-026-77458-4)、[Schroeder 2026, bioRxiv](https://www.biorxiv.org/content/10.64898/2026.04.08.717130v1)
- Doublet：[Germain 2022, F1000Research](https://f1000research.com/articles/10-979)、[Cell Genomics 2024](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00176-9)
- QC：[miQC, Hippen 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8415599/)
- 归一化与特征选择：[Ahlmann-Eltze 2023, Nat Methods](https://www.nature.com/articles/s41592-023-01814-1)、[Zappia 2025, Nat Methods](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11978513/)
- 整合：[Luecken 2022, Nat Methods](https://www.nature.com/articles/s41592-021-01336-8)、[sysVI](https://link.springer.com/article/10.1186/s12864-025-12126-3)、[scIB-E, Genome Biol 2025](https://link.springer.com/article/10.1186/s13059-025-03869-z)
- 聚类后推断：[Neufeld 2024, Biostatistics](https://academic.oup.com/biostatistics/article/25/1/270/6893953)、[ClusterDE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10418557/)
- DE：[Squair 2021](https://www.nature.com/articles/s41467-021-25960-2)、[Junttila 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9487674/)、[Nguyen 2023](https://www.nature.com/articles/s41467-023-37126-3)、[GLIMES 2025](https://link.springer.com/article/10.1186/s13059-025-03525-6)、[dreamlet 2026](https://www.nature.com/articles/s41467-026-75680-8)
- 扰动预测：[Ahlmann-Eltze 2025, Nat Methods](https://www.nature.com/articles/s41592-025-02772-6)、[Nat Biotechnol 2026](https://www.nature.com/articles/s41587-026-03307-w)
- 线粒体等位基因：[MAESTER, Miller 2022 预印本](https://www.biorxiv.org/content/10.1101/2021.03.08.434450v1.full)、[Cells 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13464780/)

## 10. 复核命令

```bash
E=/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/docs/P47_scRNAseq_best_practices_evidence_1104
cat $E/check3.log $E/check4.log   # Seurat 线粒体比例、scDblFinder 阳性对照
cat $E/github_repo_metadata_1104.tsv   # 仓库元数据（截至 2026-11-04）
# 重跑：在一个可写的临时目录中运行，check1 和 check2 会向当前目录写 CSV
/mnt/ex_8T_SSD/conda_envs/sc_tut_r/bin/Rscript $E/check4_scdblfinder_positive_control.R
```
