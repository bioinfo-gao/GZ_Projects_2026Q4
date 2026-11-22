# R for Bioinformatics Tutorial（详细版）

基础版（`../R-for-Bioinformatics-Tutorial`）的**详细版**：章节结构相同，每章内容约翻倍（共约 32 小时），
重点加强了 apply 家族、真实数据中缺失值的定位与统计、参考基因组/注释版本的混乱（GRCh37 vs GRCh38、
命名风格、liftOver、注释漂移）、ID 映射与 biomaRt，以及每一步的质量核查与稳健性分析。
每章开头有"本章路线"表；HTML 页面右侧与正文开头都有本章目录。

一套面向**有编程经验读者**的 Quarto 中文教程：
R 语言特性 → tidyverse / data.table 等数据处理生态 → ggplot2 → Bioconductor →
Bulk RNA-seq → 单细胞 RNA-seq → 肿瘤突变（oncoplot/瀑布图）与生存分析 → 工程化。

全部使用真实公开数据：

| 章节 | 数据 | 来源 |
|---|---|---|
| 4–9 | airway：气道平滑肌细胞 ± 地塞米松（Himes et al. 2014） | Bioconductor `airway` |
| 8 | pasilla：果蝇 pasilla 基因敲低（Brooks et al. 2011） | Bioconductor `pasilla` |
| 10–11 | PBMC 3k（10x Genomics） | 10x 官网 |
| 12 | Kang 2018 PBMC ± IFN-β，8 名供者 | Bioconductor `muscData` |
| 13 | Nestorowa 2016 小鼠造血干/祖细胞（Smart-seq2） | Bioconductor `scRNAseq` |
| 3、14 | TCGA 肺腺癌 MC3 体细胞突变 + TCGA-CDR 临床生存数据 | `maftools::tcgaLoad()` |
| 1、2 | TCGA-LUAD 临床表（780 列，缺失值教材） | 同上 |
| 2、7、8 | GSE52778 全部 16 个样本（GEO 表型 + NCBI GRCh38 统一重处理计数） | GEOquery / NCBI GEO |
| 5 | UCSC hg38 knownGene、hg19ToHg38 chain、EnsDb v86、Ensembl BioMart (GRCh37) | Bioconductor / UCSC / Ensembl |

## 快速开始

```bash
cd R-for-Bioinformatics-Tutorial-Detailed
./install-conda.sh         # conda/miniforge 的 R：二进制安装（推荐）
# 或：Rscript install.R      # 系统 R：从 CRAN/Bioconductor 安装（30–90 分钟）
quarto preview             # 边看边跑
quarto render              # 渲染整本书到 _book/
```

也可以在 Positron / RStudio 中逐个打开 `.qmd` 交互式运行。
中间结果缓存在 `data/`（已 gitignore）。

## 章节

| 部分 | 章节 |
|---|---|
| I R 语言与数据处理 | 01 R 语言特性 · 02 tidyverse · 03 数据处理生态（data.table/plyr/arrow/duckdb）· 04 可视化 · 05 Bioconductor |
| II Bulk RNA-seq | 06 计数与质控 · 07 DESeq2 · 08 edgeR/limma 与复杂设计 · 09 功能富集 |
| III 单细胞 | 10 质控 · 11 聚类与注释 · 12 整合与伪 bulk · 13 Bioconductor 与轨迹 |
| IV 肿瘤基因组 | 14 体细胞突变、oncoplot/瀑布图、TMB、生存分析 |
| V 工程化 | 15 可复现性与性能 |

每章都可在 Positron / RStudio 中逐块交互运行；章节之间通过 `data/` 中的缓存衔接
（例如第 7 章读取第 6 章保存的 `airway_dds.rds`），请按顺序首次运行。
