# P47 数据格式实测：10x Fixed RNA Profiling (Flex) 探针式单细胞文库

记录日期：2026-11-25　｜　记录人：Claude Code　｜　性质：内部技术笔记，非客户交付物

## 0. 一句话结论

**这批数据不是常规 3′/5′ 单细胞文库，而是 10x Fixed RNA Profiling（Flex，商品名 Chromium Fixed RNA Profiling / "Gene Expression Flex"）。8 个样本 = 2 个 physical GEX library × 每个 library 内 4 个 probe barcode，所以 `00_fastq/` 里只有 2 套 R1/R2 是正确的，不是缺文件。** 测序读到的是**合成探针 DNA**，不是转录本碱基 —— 实测 88.5% 的 R2 reads 前 50 nt 与探针序列**精确逐字匹配**。

## 1. 为什么"8 个样却只有 4 对 FASTQ"——文件数是对的

### 1.1 目录实测清点

`/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/genewiz_30-1371703644_1124/30-1371703644/00_fastq/`

截至 2026-11-25 11:33，该目录内为 **8 个 `.fastq.gz` + 8 个 `.md5` sidecar（在 `md5/` 子目录）**，共 882.345 GiB。

| 检查项 | 远端清单 | 本地实测 | 判定 |
| :--- | :---: | :---: | :---: |
| 全交付目录文件数 | 380 | 380 | 一致 |
| `00_fastq/` FASTQ 数 | 8 | 8 | 一致 |
| FASTQ md5 校验 | — | 8 OK / 0 FAILED | 通过 |
| `00_fastq/` 体积 | 882.3 GiB | 882.345 GiB | 一致 |

远端清单为 `/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/remote_listing_genewiz_20261124_0017.txt`（650 行）；其中 `00_fastq/` 段落只列出这 8 个 FASTQ。下载日志 `/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/logs/download_watch.log` 末行为 `[2026-11-25 01:47:18] DONE 下载结束。本地文件 380（远端清单 380），md5 OK 8/8，FAILED 0。自动重启 0 次。`

**⇒ 文件没有缺失。"只有 4 对"来自把常规 3′ 文库的「一样本一套 FASTQ」假设套到了 Flex 上。**

### 1.2 四类文件各是什么

每个 library 一套 4 个文件，命名 `<4个样本名拼接>_S1_L001_{I1,I2,R1,R2}_001.fastq.gz`：

| 文件 | 实测读长 | 内容 | 是否 biological read |
| :--- | :---: | :---: | :---: |
| `I1` | 10 nt | i7 **library-level** sample index | 否 |
| `I2` | 10 nt | i5 **library-level** sample index | 否 |
| `R1` | 28 nt | 16 nt GEM barcode + 12 nt UMI | 否（细胞/分子标识） |
| `R2` | 90 nt | **50 nt 合成探针** + 常量 capture/linker + 10 nt **probe barcode** | 否（是探针，不是转录本） |

实测读长由每文件前 1000 条 read 统计，长度单一：R1 全为 28，R2 全为 90，I1/I2 全为 10。

**关键区分（最容易搞混的一点）：**

- **I1/I2 = library index**，区分的是"哪个测序文库"，这里只有 2 个文库，所以 I1/I2 **不能**用来拆 8 个样本。
- **R2 末端的 probe barcode = sample index**，这才是区分 4 个样本的那一层。它不在 index read 里，而在 **biological read R2 内部**。

这就是为什么 Flex 的样本数可以超过 FASTQ 套数：**样本标签藏在 R2 里，常规按文件名拆样的思路看不到它。**

### 1.3 从原始 reads 直接验证 4 样本拆分（独立互证）

取每个 pool 的 R2 前 200,000 条 read，用 Cell Ranger 自己输出的 cell barcode 后 10 nt（即 probe barcode 序列）去反查原始 reads：

| pool | library（FASTQ 前缀） | probe barcode ID | 样本 ID | probe-BC 序列 | 200k reads 中占比 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| pool 1 | `N1-ABE-N2-ABE-N1-CBE-N2-CBE` | A-A10 | N1-ABE | GAAGGTCCCA | 19.0% |
| pool 1 | 同上 | A-B10 | N2-ABE | TGCAGTAACT | 25.7% |
| pool 1 | 同上 | A-C10 | N1-CBE | AGTATGGCTG | 27.6% |
| pool 1 | 同上 | A-D10 | N2-CBE | ACGACCCATC | 23.8% |
| pool 1 | 同上 | — | 未匹配 | — | 3.9% |
| pool 2 | `N1-Control-N2-Control-CBE-Control-ABE-Control` | A-E10 | N1-Control | TGATCATGCA | 34.8% |
| pool 2 | 同上 | A-F10 | N2-Control | GAGGGCAAAT | 34.6% |
| pool 2 | 同上 | A-G10 | CBE-Control | ACTACTCGAG | 11.8% |
| pool 2 | 同上 | A-H10 | ABE-Control | AGCTGGCTAC | 15.0% |
| pool 2 | 同上 | — | 未匹配 | — | 3.8% |

**每个 pooled FASTQ 内确实且只有 4 个 probe barcode，~96% 的 reads 能归到这 4 个样本。** 两个 pool 的 probe barcode 互不重叠（pool 1 用 A-A10…A-D10，pool 2 用 A-E10…A-H10）。复现脚本：[evidence_scripts/2_probe_barcode_demux_check.py](evidence_scripts/2_probe_barcode_demux_check.py)

注：probe barcode 在 R2 中的起始位置有少量浮动（常量区长度不完全固定），上表用"在 R2 末 20 nt 内搜索该 10-mer"的方式计数，因此 3.8–3.9% 未匹配包含测序错误与位置偏移两种情形，不代表这部分 read 无效（Cell Ranger 实测 `Valid probe barcodes` 为 97.80% / 97.91%）。

## 2. R2 读到的是探针，不是转录本——实测证据

### 2.1 前 50 nt 与探针库精确匹配

取 pool 1 的 R2 前 20,000 条 read，把每条 read 的 **bases 1–50** 拿去和 `probe_set.csv` 里 54,580 条探针序列做**精确字符串比对**：

| 比对方式 | 命中数 | 占比 |
| :--- | :---: | :---: |
| 正向精确匹配探针序列 | 17,692 | **88.5%** |
| 反向互补匹配 | 0 | 0% |
| 未精确匹配（含测序错误碱基） | 2,308 | 11.5% |

探针库中所有 54,580 条探针**长度均为 50 nt**。实测 read 结构示例（`ENSG00000204120|GIGYF2|b2702e6`）：

```
R2 (90 nt): CTCTAAGCACCATTCGGGCAAGCTATCCCTGTCATCATCTATGCTCCCAC CCCATATCCTACAACCTGCATACGCGACTT GCAGTAACTC
            └────────────── 探针序列 50 nt（逐字命中）──────────────┘ └──── 常量 capture/linker ────┘ └probe BC┘
```

**88.5% 的 read 前 50 nt 与一条已知合成探针逐字节相同，且方向全部一致** —— 这是比引用 10x 官方文档更强的本地证据：R2 的信息来源是合成 DNA 探针，不是逆转录自转录本的 cDNA。剩余 11.5% 绝大多数是带测序错误的同一批探针（Q30 RNA read 实测 96.2%），不是另一类序列。复现脚本：[evidence_scripts/1_r2_read_anatomy.py](evidence_scripts/1_r2_read_anatomy.py)

### 2.2 由此得到的硬边界

Flex 的分子机制是：一对半探针（LHS/RHS，各 25 nt）在固定细胞内与目标转录本**杂交**，杂交成功则在接合点**连接（ligation）**成一条 50 nt 完整探针，随后被捕获、扩增、测序。测序仪读到的始终是**探针自己的序列**。

**⇒ 原转录本上的任何碱基变异都不会出现在测序数据里。** 因此本批数据**在原理上**不能产出：

- 编辑位点的 ref/alt allele
- 编辑率 / heteroplasmy 比例
- 任何 base-level 的 "已编辑 / 未编辑" 细胞分型

这一点与测序深度、是否索取 BAM 无关：即便 GENEWIZ 补发 BAM（`config.csv` 写了 `create-bam,true`，但本地镜像的 380 个文件中**不含任何 `.bam`/`.cram`**），BAM 里比对的也是探针序列。**索取 BAM 不能解除这条限制。**

## 3. 探针面板对线粒体基因组的覆盖——比"不能测碱基"更严的限制

面板为 `Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv`：54,580 条探针 / 19,074 个基因；其中 53,444 条 `included=TRUE`，1,136 条 `FALSE`（被排除出 filtered probe set）。

### 3.1 chrM 只有 12 个基因、每基因只有 1 条探针

| 项目 | 实测 |
| :--- | :---: |
| 面板内 `MT-` 前缀基因数 | **12**（全部为蛋白编码） |
| mt-tRNA（`MT-T*`）探针数 | **0** |
| mt-rRNA（`MT-RNR1/2`）探针数 | **0** |
| 每个 MT 基因的探针数 | **1**（全部 `region=unspliced`） |
| 其中 `included=FALSE` | MT-CO1（1 个） |

对比全面板：15,309 个基因有 3 条探针、1,011 个有 2 条、1,830 个有 1 条。**MT 基因属于"单探针"那一类，这对后续判读有直接后果（见 3.3）。**

### 3.2 12 条 MT 探针在 chrM 上的足迹

把 12 条探针序列比对到 `GRCh38.primary_assembly.genome.fa` 的 chrM（实测长度 16,569 bp，与 rCRS 一致），**12 条全部唯一定位**：

| 基因 | 是否 included | 链 | chrM 足迹（1-based） | 连接点（ligation junction） |
| :--- | :---: | :---: | :---: | :---: |
| MT-ND1 | TRUE | − | 3,309–3,358 | 3,333 / 3,334 |
| MT-ND2 | TRUE | − | 5,122–5,171 | 5,146 / 5,147 |
| MT-CO1 | **FALSE** | − | 7,121–7,170 | 7,145 / 7,146 |
| MT-CO2 | TRUE | − | 7,866–7,915 | 7,890 / 7,891 |
| MT-ATP6 | TRUE | − | 9,050–9,099 | 9,074 / 9,075 |
| MT-CO3 | TRUE | − | 9,743–9,792 | 9,767 / 9,768 |
| MT-ND3 | TRUE | − | 10,305–10,354 | 10,329 / 10,330 |
| MT-ND4L | TRUE | − | 10,645–10,694 | 10,669 / 10,670 |
| MT-ND4 | TRUE | − | 10,897–10,946 | 10,921 / 10,922 |
| MT-ND5 | TRUE | − | 12,429–12,478 | 12,453 / 12,454 |
| MT-ND6 | TRUE | + | 14,361–14,410 | 14,385 / 14,386 |
| MT-CYB | TRUE | − | 15,713–15,762 | 15,737 / 15,738 |

机器可读 SSOT（含探针序列）：[chrM_probe_footprints.csv](chrM_probe_footprints.csv)，12 行。复现脚本：[evidence_scripts/3_mt_probe_footprint.py](evidence_scripts/3_mt_probe_footprint.py)、[evidence_scripts/4_write_mt_probe_footprint_csv.py](evidence_scripts/4_write_mt_probe_footprint_csv.py)

**覆盖总量：12 × 50 = 600 bp，占 chrM 16,569 bp 的 3.62%；只算 `included=TRUE` 的 11 条为 550 bp = 3.32%。**

### 3.3 一个仅存的间接信号，以及它为什么在本批数据里基本不可用

**机制上存在**一条间接路线：若编辑后的碱基落在探针足迹内（尤其跨过连接点），杂交/连接效率下降 → 该探针的 UMI 下降 → 表现为该基因"表达降低"。Cell Ranger 甚至为此报告 `Reads half-mapped to probe set`（连接失败的读段），本批 8 个样本实测在 0.147%–0.169% 之间。

**但在本批数据里这条路线有三重阻断：**

| 阻断 | 实测依据 | 后果 |
| :--- | :---: | :---: |
| ① 落点概率低 | 探针只覆盖 chrM 的 **3.62%** | 编辑位点先验上有 ~96% 概率根本不在任何探针下 |
| ② 经典靶点全部未覆盖 | m.3243 / m.3460 / m.8344 / m.8993 / m.11778 / m.14484 **无一落在 12 个足迹内** | 若靶点是常见致病位点（多为 mt-tRNA），面板**零覆盖** |
| ③ 无基因内对照 | MT 基因**每基因仅 1 条探针** | "探针信号下降"与"该基因表达真的下降"**完全共线**，无法区分 |

第 ③ 条是方法学上最要命的：对有 3 条探针的核基因，可以用"同基因其他探针作内参"的 probe-ratio 把"杂交失败"从"表达变化"里分离出来；**MT 基因没有这个内参**。所以即使编辑位点恰好落在足迹内，单探针下降也**不能**被定性为编辑证据。

**⇒ 结论：这条间接路线只能作为 hypothesis-generating 的探索性观察，在任何情况下都不构成编辑成功的证据。而且必须先确认靶点坐标是否落在 `chrM_probe_footprints.csv` 的 600 bp 内；若不在，这条路线连探索都无从开始。**

**客户给出靶点坐标后的 10 秒判定方法**（可直接复制运行）：

```bash
conda run -n regular_bioinfo python3 -c "
import csv
POS = 3243   # ← 改成客户提供的 chrM 1-based 坐标
rows=list(csv.DictReader(open('/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/data_format/chrM_probe_footprints.csv')))
hit=[r for r in rows if int(r['chrM_start_1based'])<=POS<=int(r['chrM_end_1based'])]
print(f'chrM:{POS} ->', hit[0]['gene_name']+' '+hit[0]['ligation_junction_between'] if hit else 'NOT covered by any probe (3.62% 覆盖率外)')
"
```

## 4. Flex 特有的两个分析期陷阱（与常规 scRNA 不同）

### 4.1 cell barcode 是 26 nt，不是 16 nt

Cell Ranger 输出的 barcode 形如 `AAACCAATCCAATCATGAAGGTCCCA-1`，结构为 **16 nt GEM barcode + 10 nt probe barcode**。**同一个 GEM 液滴里可以同时存在多个样本的细胞**（Flex 是固定后混样再上机）。后果：

`multiplexing_analysis/frp_gem_barcode_overlap.csv`（各 7 行）实测，同 pool 内任意两个 probe barcode 的 GEM 重叠率为 **11.4%–17.4%**：

| pool | 重叠率范围 | 最低 / 最高 |
| :--- | :---: | :---: |
| pool 1（edited） | 13.5%–17.1% | A-A10×A-D10 / A-C10×A-D10 |
| pool 2（control） | 11.4%–17.4% | A-G10×A-H10 / A-E10×A-H10 |

**⇒ 同 pool 内存在样本间的 ambient / multiplet 交叉污染，且是结构性的（不是异常）。做组间比较时，同 pool 的样本之间并非完全独立观测；ambient 背景校正必须按 pool 内多样本情形处理，不能当成单样本的 soup。**

### 4.2 不能把常规 STARsolo 流程套用于 Flex

项目 `scripts/` 下已有的 `2_build_starsolo_ref.sh`、`3_detect_chemistry.py`、`4_starsolo.sh` 是针对常规 10x v2/v3/v4 UMI 文库写的，**对 Flex 不适用**，原因是实测结构性的、不是参数能调的：

| 常规流程的假设 | Flex 的实测事实 |
| :--- | :---: |
| cell barcode 16 nt，全部在 R1 | 26 nt，其中 10 nt 在 **R2** 内部 |
| R2 是转录本 cDNA，需比对到基因组/转录组 | R2 是 50 nt 合成探针，需比对到 **probe set** |
| 用 10x whitelist 校正 barcode | 需同时校正 GEM whitelist **和** 384 个 probe barcode |
| 计数单位 = 基因 | 计数单位 = **探针**，再按基因汇总 |

GENEWIZ 已用 Cell Ranger 8.0.0 + 对应 probe set 完成拆样与定量，产物齐全（见 §5）。**在没有经验证的 FRP 重计数方案前不重跑 STARsolo；这些脚本应标注为"仅适用于常规 3′/5′ 文库，本项目不使用"，而不是删除。**

## 5. 公司已交付产物清单与质量指标

### 5.1 产物结构

`01_analysis/` 下 2 个 library 目录，各含 library 级产物 + `per_sample_outs/<sample_id>/` 4 个样本目录。**单样本级产物齐全**：

| 产物 | 路径模式 | 说明 |
| :--- | :---: | :---: |
| 样本 filtered 矩阵 | `per_sample_outs/<s>/sample_filtered_feature_bc_matrix{,.h5}` | 公司 cell-called 集合，分析起点 |
| 样本 raw 矩阵 | `per_sample_outs/<s>/sample_raw_feature_bc_matrix{,.h5}` | cell-calling / ambient 背景核查用 |
| **探针级矩阵** | `per_sample_outs/<s>/sample_raw_probe_bc_matrix.h5` | **54,580 探针 × barcode**，probe-ratio 分析的唯一入口 |
| 分子信息 | `per_sample_outs/<s>/sample_molecule_info.h5` | 饱和度 / 重抽样 |
| 样本指标 | `per_sample_outs/<s>/metrics_summary.csv` | 单样本 QC |
| 公司聚类/DE | `per_sample_outs/<s>/analysis/{clustering,diffexp,pca,tsne,umap}` | **cluster-wise，不是样本级推断** |
| 面板定义 | `<library>/probe_set.csv` | 54,580 行 |
| library 指标 | `<library>/qc_library_metrics.csv` | pool 级 |
| 多重化诊断 | `<library>/multiplexing_analysis/{cells_per_tag.json,frp_gem_barcode_overlap.csv}` | §1.3、§4.1 数据来源 |

实测 `sample_raw_probe_bc_matrix.h5`（N1-ABE）：54,580 探针 × 732,922 barcode，258,774,469 个非零值；`features/` 内含 `probe_region`、`filtered_probes`、`target_sets` 字段。

`cells_per_tag.json` 列全 384 个 probe barcode 位，每个 library 内**仅 4 个非空**，细胞数与 QC 表逐一吻合。

### 5.2 library 级指标

| 指标 | pool 1（edited） | pool 2（control） |
| :--- | :---: | :---: |
| Number of reads | 7,269,451,607 | 6,478,083,348 |
| Cells（library 合计） | 70,704 | 65,340 |
| Sequencing saturation | 33.3% | 26.7% |
| Valid barcodes | 96.43% | 96.54% |
| Valid probe barcodes | 97.80% | 97.91% |
| Confidently mapped reads in cells | 97.77% | 97.69% |
| Q30 RNA read | 96.18% | 96.23% |
| Estimated UMIs from genomic DNA | 0.0088% | — |

### 5.3 样本级指标（公司值，非本项目独立重算）

| 样本 ID | pool | probe BC | Cells | Median UMI/cell | Median genes/cell | Reads in cells | Half-mapped |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| N1-ABE | 1 | A-A10 | 15,245 | 47,228 | 8,625 | 1,358,931,130 | 0.164% |
| N2-ABE | 1 | A-B10 | 19,027 | 52,216 | 8,867 | 1,859,024,169 | 0.147% |
| N1-CBE | 1 | A-C10 | 20,706 | 48,541 | 8,878 | 1,941,921,733 | 0.151% |
| N2-CBE | 1 | A-D10 | 15,726 | 58,424 | 9,175 | 1,693,758,770 | 0.169% |
| N1-Control | 2 | A-E10 | 20,780 | 61,946 | 9,109 | 2,222,181,629 | — |
| N2-Control | 2 | A-F10 | 20,053 | 61,535 | 9,199 | 2,192,306,043 | — |
| CBE-Control | 2 | A-G10 | 12,480 | **26,250** | **6,840** | 746,751,973 | — |
| ABE-Control | 2 | A-H10 | 12,027 | 42,071 | 8,209 | 949,000,266 | — |

合计 136,044 cells（12,027–20,780/样本）。中位 6,840–9,199 genes/cell，对 293T 属高质量。

**两处必须在分析时处理的不均衡：**

1. **每样本 read 分配极不均衡**：`Reads in cells` 从 7.47e8（CBE-Control）到 2.22e9（N1-Control），**差 2.98 倍**。这是 probe barcode 间上样量/捕获效率差异造成的，不是生物学差异。⇒ 样本间比较 UMI/genes per cell 之前必须先做深度校正或下采样，**不能直接把 CBE-Control 的低复杂度读成生物学效应**。
2. **CBE-Control 是复杂度最低的样本**（26,250 UMI / 6,840 genes，均为 8 样本最低），而它是**对照组**成员。⇒ 若后续出现"CBE 相关对比差异最大"的结果，必须先排除这是深度/质量 artifact。

## 6. 本文件的事实来源与复核方式

| 本文主张 | 来源 | 复核方式 |
| :--- | :---: | :---: |
| 文件数 380/380、md5 8/8 | 下载日志 + 远端清单 | `tail -1 /mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/logs/download_watch.log` |
| 读长 28/90/10/10 | 每文件前 1000 read 实测 | `zcat <f> | head -4000 | awk 'NR%4==2{print length($0)}' | sort | uniq -c` |
| R2 前 50 nt 88.5% 精确命中探针 | 本地实测 | `conda run -n regular_bioinfo python3 -I evidence_scripts/1_r2_read_anatomy.py` |
| 每 pool 内 4 个 probe barcode、~96% 归属 | 本地实测 | `conda run -n regular_bioinfo python3 -I evidence_scripts/2_probe_barcode_demux_check.py` |
| 12 条 MT 探针的 chrM 足迹、覆盖 3.62% | 本地比对 GRCh38 chrM | `conda run -n regular_bioinfo python3 -I evidence_scripts/3_mt_probe_footprint.py` |
| GEM 重叠 11.4–17.4% | 公司诊断文件 | `cat 01_analysis/*/multiplexing_analysis/frp_gem_barcode_overlap.csv` |
| 样本级 QC 数值 | 公司 `metrics_summary.csv` | `grep -E ",Cells,|Median " 01_analysis/*/per_sample_outs/*/metrics_summary.csv` |

**我的判断与推理（非实测，需用户确认或在报告中标为推断）：**

- "R2 剩余 11.5% 未精确匹配主要是测序错误" —— 依据是 Q30 RNA read 96.2% 与 `Reads mapped to probe set` 98.4%，但未逐条做容错比对确认。
- "每样本 read 分配不均衡源于上样量/捕获效率差异" —— 机制推断；湿实验投入量未知。
- §3.3 的 probe-ratio 内参思路对本面板 MT 基因不可用是**结构性事实**（1 探针/基因，实测）；但"落在连接点上影响更大"是机制推断，本批数据未验证。
