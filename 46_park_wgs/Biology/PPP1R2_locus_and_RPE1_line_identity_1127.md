# 靶点与材料的身份 —— PPP1R2 基因/蛋白、这个位点的生物信息学特征，以及三层细胞系命名

- **创建日期（不可更改）**：2026-11-27
- **物种 / 参考**：*Homo sapiens*，GRCh38，注释用 `gencode.v45.annotation.gtf`
- **涉及项目**：P31（2 样本）、P46（10 样本）—— 两轮的靶点逐字相同
- **写这一份的起因**：PPP1R2 在 P31/P46 的文档里一直只以**坐标**出现（`chr3:195,542,947 / 195,543,149`），
  从没写过**这是个什么基因、删掉的 202 bp 在基因结构里是什么位置、删掉之后这条等位还剩什么**。
  本文补这一层，并顺带把 `sLK_of_CRISPR_Editing/dual_cut_per_cut_genotyping_1126.md` §5
  明确挂着未做的那个核对（Sample_5 第二条等位是否移码）做完。
- **2026-11-27 增补 §2A**：材料命名。`hTERT RPE-1` / `RPE-1_53bp1_mcherry` / `RPE1_WT_Ctr`
  这三个名字此前散落在送样单、计划文档与 P31 文档里，**没有一处说明它们是同一条谱系的三层**，
  而其中 `WT` 与 `Ctr` 两个词都名不副实。§2A 把三层关系、每个 token 的含义与四个误读点写在一处。
- **定位**：本文回答的是"**这些名词到底指什么**" —— 靶点基因、靶点位点、以及材料与样本名。
  内容是**注释层面与命名层面的事实**，**不是实验功能验证**：
  凡涉及"蛋白还有没有功能"，一律只到注释预测为止。

## ⚠ 来源标注约定（全文强制）

| 标记 | 含义 |
| :--- | :---: |
| **【客】** | 客户书面陈述 |
| **【测】** | 我方实测，本文给出命令或产物路径 |
| **【推】** | 我方推断 |
| **【景】** | 领域背景，**不作为本项目证据使用** |
| **【景·待核】** | 凭印象的领域知识，**写进任何对外文字之前必须查原始文献** |

---

## 1. 三段话说完

**第一段 —— 这个基因是什么。**【景】+【测】
PPP1R2 编码 **protein phosphatase 1 (PP1) inhibitor-2**，习称 **I-2**（旧名 IPP-2）：
一个 **205 aa**【测】、内在无序（intrinsically disordered）、耐热耐酸的小蛋白，
功能是**结合并抑制 PP1 催化亚基**，同时参与 PP1 的折叠成熟与分裂期调控。
它**不是**客户的研究对象 —— 客户要的是在一个指定位置制造双链断裂来观察 53BP1，
PPP1R2 只是"那个位置"。**所以本基因的生物学在本项目里的作用是潜在混杂因素，不是读出指标。**

**第二段 —— 删掉的 202 bp 恰好是这个基因的翻译起点，而且删得异常"整齐"。**【测】
缺失 `chr3:195,542,948–195,543,149` 落在 **exon 1**，删去
**the 3′-most 124 bp of the 5′UTR（紧接 ATG 上游那一段）+ 起始密码子 ATG + 密码子 1–26（78 bp CDS）**，
124 + 78 = 202 ✓。
**78 bp 正好是 26 个完整密码子** ⇒ 残余下游 CDS **仍在原读框上**，断点正落在密码子 26/27 之间。
**转录起点（TSS）与启动子未被删** —— 保留了 176 bp 的 5′UTR，所以这不是启动子缺失，
**是"起始密码子 + N 端 26 残基"缺失**，转录大概率照常进行。

**第三段 —— 两条等位的后果都算清楚了，但都只到注释层面。**【测】+【推】

| 等位 | 分子事件 | 注释层面的预测后果 |
| :--- | :---: | :---: |
| **Δ202**（Sample_1–7 携带的那条） | ***PPP1R2*** **exon 1**，`chr3:195,542,948–195,543,149`（202 bp）：删除起始密码子 ATG + 密码子 1–26 | 原起始密码子消失；扫描模型下第一个可用 AUG 开出 **156 aa 的 N 端截短产物（残基 50–205）**，前面还垫着两个 5 aa / 2 aa 的小 ORF |
| **Sample_5 的第二条** | ***PPP1R2*** **exon 1**，`chr3:195,542,947–195,542,954` **8 bp 缺失**（CIGAR 实测） | **移码**：残基 1–24 正常，自第 25 位起错读，**45 aa 处提前终止**（野生型 205 aa） |

⇒ **Sample_5 的"双等位"在分子层面确实是两条都被改动**，且两条都预测为功能丧失 ——
这一条此前被标为"未核对、不下结论"，**现在可以下结论了**（§5.3）。


---

## 2. 缩写表

| 缩写 | full name | 是什么 |
| :--- | :---: | :---: |
| **ATCC** | American Type Culture Collection | 细胞系保藏机构；hTERT RPE-1 编号 CRL-4000 |
| **CDS** | coding sequence | 编码序列 |
| **CIGAR** | Compact Idiosyncratic Gapped Alignment Report | SAM/BAM 里记录比对如何分段（M/I/D/S 等）的那一列 |
| **DSB** | double-strand break | DNA 双链断裂 |
| **EGFP** | enhanced green fluorescent protein | 增强型绿色荧光蛋白 |
| **gRNA** | guide RNA | 向导 RNA |
| **GSK-3** | glycogen synthase kinase 3 | 糖原合酶激酶 3 |
| **I-2** | inhibitor-2 | PPP1R2 的蛋白产物习用名 |
| **IDP** | intrinsically disordered protein | 内在无序蛋白 |
| **MANE** | Matched Annotation from NCBI and EMBL-EBI | NCBI 与 Ensembl 共同指定的代表转录本 |
| **MAPQ** | mapping quality | 比对质量值；越低表示该读段的定位越不唯一 |
| **NHEJ** | non-homologous end joining | 非同源末端连接 |
| **NMD** | nonsense-mediated decay | 无义介导的 mRNA 降解 |
| **PP1** | protein phosphatase 1 | 丝氨酸/苏氨酸蛋白磷酸酶 1；`PP1c` 指其催化亚基 |
| **PTC** | premature termination codon | 提前终止密码子 |
| **RNP** | ribonucleoprotein | 核糖核蛋白；此处指 Cas9 蛋白与 gRNA 的体外预组装复合物 |
| **TSS** | transcription start site | 转录起点 |
| **uORF** | upstream open reading frame | 主 ORF 上游的小开放读框 |
| **UTR** | untranslated region | 非翻译区 |

---

## 2A. 材料命名：三个名字指的是同一条谱系的三层

> **这一节回答的问题**：`hTERT RPE-1`、`RPE-1_53bp1_mcherry`、`RPE1_WT_Ctr` 各是什么意思。
> **一句话答案**：**它们不是三株细胞，是同一条谱系的三层**；而且后两个名字里的
> `WT` 与 `Ctr` **都名不副实**，是本项目最容易被外人误读的两个词。

### 2A.1 三层谱系

| 层级 | 名字（原始写法与出处） | 是什么 | 谁用这个名字 |
| :--- | :---: | :---: | :---: |
| **① 亲本细胞系** | **`hTERT RPE-1`**（ATCC **CRL-4000**）【客·问卷】 | 人视网膜色素上皮细胞系，靠导入 hTERT 永生化（**非病毒癌基因**）；近二倍体、p53 通路完整 | P31 与 P46 的**共同亲本** |
| **② 报告器系** | 送样单 `cell_line` 字段：**`RPE-1_53BP1-FFR-mCherry monoclonal`**【客】；计划文档另写作 `RPE-1_53BP1-mCherry`；口语写法 `RPE-1_53bp1_mcherry` | 在 ① 上用**逆转录病毒**整合报告器 `mCherry-BP1-2 pLPC-Puro`（全长 **7,985 bp**【测】），荧光分选后 3 cells/well 扩增 | P46 送样单 **10/10 样本**的 `cell_line` 字段 |
| **③ 单克隆（= 实际样本）** | P31：**`RPE1_WT_Ctr`**、`RPE1_cg13_Single`【客】；P46：`Ctr_1–3`、`Sample_1–7`【客】 | 在 ② 上做 caged dual-cut CRISPR 之后挑出的**单细胞来源克隆**，传代均 < 20 | 各轮的样本名 |

⇒ **③ 是送来测序的东西，② 是 ③ 的亲本，① 是 ② 的亲本。**
⇒ `RPE1_WT_Ctr` 与 `Ctr_1–3` 处在**同一层**（都是第 ③ 层的对照克隆），
只是分属两轮：**`RPE1_WT_Ctr` 是 P31 的样本名，P46 不用它**。

### 2A.2 名字逐个 token 展开

| token | full name | 含义 |
| :--- | :---: | :---: |
| **hTERT** | human telomerase reverse transcriptase | 人端粒酶逆转录酶 —— **永生化手段**，不是基因型描述 |
| **RPE** | retinal pigment epithelium | 视网膜色素上皮，组织来源 |
| **-1** | — | **亚系编号**，不是"1 个拷贝" |
| **53BP1** | p53-binding protein 1 | DSB 修复因子（基因名 `TP53BP1`，chr15） |
| **FFR** | focus-forming region | 成灶区，53BP1 能在断裂处聚成 focus 的最小片段；构建体内为 **aa 1225–1716**【测】 |
| **BP1-2** | — | 构建体内对该 53BP1 片段的编号，载体名 `mCherry-BP1-2` 的后半 |
| **mCherry** | — | 红色荧光蛋白；**我方载体定量实际依赖的就是这一段**（FFR 那半因与内源 `TP53BP1` 同源被屏蔽） |
| **pLPC** | — | 逆转录病毒载体骨架 |
| **Puro / PuroR** | puromycin resistance | 嘌呤霉素抗性筛选标记 |
| **WT** | wild type | ⚠ 见 2A.3 ① —— **仅指靶点，不指基因组** |
| **Ctr** | control | ⚠ 见 2A.3 ② —— **不是 CRISPR-naive** |
| **cg13** | — | 客户对 P31 那个被编辑克隆的编号 |
| **Single** | single-allelic | P31 克隆的**先验**基因型（凝胶/Sanger），非测序结论 |

⇒ **拼起来读**：`RPE-1_53BP1-FFR-mCherry` = "带 53BP1 成灶区–mCherry 融合报告器的 RPE-1"；
`RPE1_WT_Ctr` = "RPE-1 背景、靶点未见缺失、充当对照的那个克隆"。

### 2A.3 ⚠ 四个会被误读的地方

**① `WT` 只对靶点成立，对基因组不成立。**【客】+【测】
`RPE1_WT_Ctr` 携带**两个**整合构建体 —— 客户原文自述其
"contains both Cas9-EGFP and 53BP1-FFR-mCherry"。
⇒ 它对 *PPP1R2* 是 WT/WT，但**绝不是"野生型细胞"**。
写进任何文字时，`WT` 后面必须跟上限定语（"WT **at the target locus**"）。

**② `Ctr` 不是 CRISPR-naive —— 这是客户自己先说的。**【客】
客户原话：
> **So this one is not ideal as a control** but gDNA from this cell line was the only monoclonal one
> without single-allelic deletion.

P46 的 `Ctr_1–3` 处境完全相同：同一载体、同样经过 caged dual-cut CRISPR 处理，
**只是凝胶上未见缺失**。
⇒ **这正是全项目结论层级止于"与编辑共分离（segregating with the edit）"而非"由编辑导致"的根源之一。**
（另一个根源是两组核型不同源，见 `hTERT_RPE1_cell_line_biology_1125.md` §4.5。）

**③ 第 ② 层的名字并不决定 Cas9 状态 —— 同名材料实测不同。**【测】

| 轮次 | 样本 | Cas9 CDS 实测拷贝数 | 含义 |
| :--- | :---: | :---: | :---: |
| **P31** | `RPE1_WT_Ctr` | **1.17** | 带**整合的** Cas9-EGFP |
| **P31** | `RPE1_cg13_Single` | **0.00** | 不带 |
| **P46** | 全部 10 个克隆 | **0.000–0.235** | 均不带（本轮 Cas9 走 **RNP**，未转质粒【客】） |

⇒ **两轮的 `cell_line` 字段写的是同一个名字，Cas9 整合状态却不同。**
**名字不能当基因型用** —— 载体状态每轮必须实测。
P46 数字来源：`analysis/stage3_vector_v2_20261124/tables/vector_coverage_n10.tsv`
（`region == "Cas9"`，取 `est_copies_per_diploid` 列）；P31 数字见
`P31_P46_caged_dual_cut_wgs_rationale_1125.md` §6.5。

**④ 写法在文档之间漂移 ⇒ 检索时必须放宽。**【测】
同一个东西至少有这些写法：`RPE-1` / `RPE1`（连字符有无）、`53BP1` / `53bp1`（大小写）、
`-` / `_`（分隔符）、以及**`FFR` 被省掉**的 `RPE-1_53BP1-mCherry`。
⇒ 在项目里找这些名字，**必须 `grep -i` 且连字符与下划线两种都试**，
否则会得到"没有这个东西"的假阴性。

**样本名的 SSOT 是 `scripts/0_sample_table.csv`**：
用 `safe_name` 列（`Ctr_1`…`Sample_7`）做一切脚本与路径；
送样单的原始写法（`Sample 1` 带空格、`Ctr_2 ` 带尾随空格）保留在同表的
`sample_name` 与 `notes` 列。**本文不复制这张表** —— 以该 CSV 为准。

### 2A.4 对报告的含义

交付报告 1127 版 header 写的是 **`hTERT RPE-1 cell line, single-cell–derived clones`**，
即**第 ① 层 + 克隆来源**；报告器则在 Methods 与载体一节另行交代。
⇒ **这不是错**（客户问卷确实写了 hTERT RPE-1 / ATCC），
但只读 header 的人**不会知道这些克隆带着一个整合的报告器**。
是否把 header 的 Tissue / Cell 改为送样单 `cell_line` 字段的原始写法，
列入 §9 待办 **#6** 由用户拍板 —— 报告正处内部审阅阶段，本文不擅自改。

---

## 3. 基因本身（全部【测】，来自 GENCODE v45）

| 项 | 内容 |
| :--- | :---: |
| 基因 ID | `ENSG00000184203.8`，HGNC:9288 |
| 位置 / 链 | **chr3:195,514,428–195,543,386**，**负链**（⇒ 坐标越大越靠 5′） |
| 跨度 | 28,959 bp |
| 转录本数 | **7**（3 protein_coding + 2 NMD + 2 retained_intron；其中 `ENST00000438848` 带 `cds_end_NF`，3′端不完整） |
| MANE Select | **`ENST00000618156.5`（PPP1R2-206）**，6 exon，CDS 615 bp，**205 aa**，CCDS3309 |
| 寄宿 ncRNA | `RNU6ATAC24P`（snRNA，chr3:195,523,087–195,523,221）落在内含子内 ⇒ 基因带 `ncRNA_host` 标签 |

**MANE 转录本的 exon 1 结构**（负链，坐标降序 = 转录方向）：

| 要素 | 坐标 | 长度 | 是否落在本项目缺失内 |
| :--- | :---: | :---: | :---: |
| exon 1 5′端（TSS） | 195,543,325 | — | **否**（保留） |
| 5′UTR | 195,543,026–195,543,325 | 300 bp | **部分**：被删的是 **the 3′-most 124 bp of the 5′UTR**（`195,543,026–195,543,149`，紧接 ATG 上游）；保留的是 **the 5′-most 176 bp**（`195,543,150–195,543,325`，含 TSS）

原文CC 的错误写法： 5′UTR	195,543,026–195,543,325	300 bp	部分：远端(distal)124 bp被删，近端(proxima) 176 bp 保留， CC的意思是离TSS远，这不符合转录本表述惯例 |
| **起始密码子 ATG** | **195,543,023–195,543,025** | 3 bp | **是** |
| exon 1 内 CDS | 195,542,904–195,543,025 | 122 bp | **部分**：78 bp 被删，44 bp 保留 |
| exon 1 3′端 | 195,542,904 | — | 否 |

⚠ **本文一律不用 distal / proximal（更不用"远端 / 近端"）描述这个 5′UTR。**
**这两个词本身不含锚点**，而本例有三个互相冲突的候选锚点：

| 锚点 | 来自哪种惯例 | 按它读，"distal 124 bp" 指哪一段 |
| :--- | :---: | :---: |
| **TSS / promoter** | promoter 文献（`proximal promoter` / `distal enhancer`） | `195,543,026–195,543,149`（**离 TSS 远、紧接 ATG**）✓ 这是本文原意 |
| **ATG / CDS** | transcript 结构描述，多数读者的默认 | `195,543,150–195,543,325`（**离 ATG 远、含 TSS**）✗ 与原意相反 |
| **基因组坐标** | 看 IGV 时的直觉 | ***PPP1R2*** **在负链** ⇒ 坐标方向与转录方向相反，再反一次 |

⇒ 全文改用 **transcript 解剖学方位**：`the 5′-most … of the 5′UTR` / `the 3′-most … of the 5′UTR`，
后面附坐标区间。**这类表述不需要锚点**，在正负链上读法一致。
本例：被删 = `the 3′-most 124 bp`（紧接 ATG 上游）；保留 = `the 5′-most 176 bp`（含 TSS）。

⚠ **起始密码子位置在全部 4 个有注释起始密码子的转录本中完全一致**
（`ENST00000618156` / `438848` / `625807` / `413183`，均为 195,543,023–195,543,025）
⇒ **这个缺失一次命中全部可翻译的转录本**，不存在"换个转录本就躲过去"的情形【测】。
（第 5 个带 CDS 的转录本 `ENST00000418939` 标为 `cds_start_NF`，即注释本身未定出其起始密码子，
不构成一条独立的逃逸路径。）

---

## 4. 蛋白的生物学：I-2 是什么角色

> 本节全部为领域背景，**不作为本项目证据**。标【景·待核】的几条凭印象，对外引用前必须查原始文献。

### 4.1 核心功能【景】

PP1 是真核细胞最主要的 Ser/Thr 磷酸酶之一，但它的催化亚基 PP1c **几乎从不单独行动** ——
它靠几十种调节亚基决定去哪里、对谁去磷酸化。**I-2 是其中最早被发现的抑制亚基之一**：
结合 PP1c 形成无活性复合物。

**两个层面的作用**：
1. **抑制**：占据 PP1c 的底物结合沟与活性位点附近，使其失活；
2. **伴侣 / 成熟（chaperone）**：新合成的 PP1c 需要被扣住、正确折入金属中心，I-2 参与这个过程。
   ⇒ 这解释了为什么 I-2 的缺失不等于"PP1 活性简单上升"，**后果可以是双向的**。

### 4.2 结构与关键基序【景·待核】

| 要素 | 内容 | 在本项目缺失片段内？ |
| :--- | :---: | :---: |
| **`KGILK` 基序** | 文献所述 I-2 的 N 端 PP1 对接基序（占据 PP1c 的 RVxF 疏水沟） | **是** —— 实测位于残基 **12–16**【测】 |
| **`PSTP` 中的 Thr** | GSK-3 磷酸化位点，经典文献写作 Thr72（本参考序列中为残基 **73**，`P71-S72-T73-P74`【测】） | 否（保留） |
| 整体构象 | 内在无序，SDS-PAGE 表观分子量高于理论值 | — |

⚠ **编号须注意**：文献中的 "Thr72" 来自早期（兔/人）序列编号，**与本参考序列的 205 aa 编号差 1 位**。
任何写进对外文字的残基编号，必须注明所用参考转录本。

### 4.3 分裂期相关【景·待核】

I-2 被报道定位于中心体与纺锤体，参与染色体分离，并有文献报道其**激活 Aurora A**。
⇒ 若属实，**在一个靠 DNA 损伤应答计时的成像体系里，敲掉 I-2 不是完全中性的背景改动**（见 §6.2）。
**这条目前只是印象，是本文最需要查文献的一条。**

---

## 5. 本项目靶点的注释层面后果

### 5.1 缺失的精确构成（【测】，可复算见 §8）

```
负链基因，转录方向 = 坐标降序
                                   删除区 202 bp
          |<-- 保留 176 bp 5'UTR -->|<------------------------>|<- 保留 44 bp CDS ->
  TSS                              CUT2                      CUT1
195,543,325 ................. 195,543,149 ......... 195,542,948|195,542,947 ... 195,542,904
                                            ATG
                                     195,543,023-025
                                   (密码子 1-26 在其下游)
```

| 被删掉的 | 长度 |
| :--- | :---: |
| the 3′-most part of the 5′UTR（紧接 ATG 上游，`195,543,026–195,543,149`） | 124 bp |
| 起始密码子 + 密码子 2–26 | 78 bp（= **26 个完整密码子**） |
| **合计** | **202 bp** ✓ |

**被删去的 26 个残基**【测】：`MAASTASHRPIKGILKNKTSTTSSMV`
⇒ 含 §4.2 的 `KGILK`（残基 12–16）。

**两个切点在基因结构中的位置不对称**，这一点对解读有用【测】：

| 切点 | 坐标 | 落在 | 后果 |
| :--- | :---: | :---: | :---: |
| **CUT1** | 195,542,947 | **exon 1 的 CDS 内**，密码子 26/27 之间 | 此处的任何 indel 都直接改读框 |
| **CUT2** | 195,543,149 | **5′UTR 内** | 此处的小 indel 不改读框，只动 UTR 长度 |

⇒ **只有 CUT1 处的 indel 会产生编码后果。** Sample_5 的第二条等位正是 CUT1 处事件（§5.3）。

### 5.2 Δ202 等位：还能不能翻译出东西（扫描模型预测，【测】+【推】）

把 Δ202 转录本重建出来（保留的 176 bp 5′UTR + 44 bp exon-1 残段 + exon 2–6），逐个 AUG 扫描：

| 转录本内位置 | ORF 长度 | 产物 |
| :--- | :---: | :---: |
| nt 202 | 5 aa | `MSTRS*` —— **uORF** |
| nt 241 | 2 aa | `MK*` —— uORF |
| **nt 246** | **156 aa** | **`MNILATYHPAD...` = 野生型残基 50–205**（N 端截短，**缺 KGILK**） |
| nt 294 | 140 aa | 野生型残基 66–205 |

**可以说的**【推】：Δ202 删掉了正规起始密码子；若核糖体扫描越过两个小 uORF 后在 nt 246 起始，
理论上可产出一个**缺少 N 端 49 个残基、因而缺少 PP1 对接基序的截短蛋白**。
**不能说的**：这是否真的发生 —— 需要 RNA / 蛋白数据，本轮是 WGS，**原理上回答不了**。
⇒ 对客户的措辞只能到："**该等位的正规翻译起点被删除**；是否存在替代起始，本轮数据无法判定。"

### 5.3 ✅ 补做：Sample_5 第二条等位 —— 移码，45 aa 截短

> 这一条回答 `dual_cut_per_cut_genotyping_1126.md` §5 结尾标注的"本文未做这个核对"。

**缺失坐标由 CIGAR 直接定出**【测】：读段 `LH00453:833:25FGWWLT4:5:1460:14667:29402`
为 `138M8D12M` @ POS 195,542,809 ⇒ 删除 ***PPP1R2*** **exon 1 内的 `chr3:195,542,947–195,542,954`，8 bp**
（另一条读段 `135M8D15M` @ 195,542,812 给出同一区间，互证）。
该区间完全落在 exon 1 的 CDS（`chr3:195,542,904–195,543,025`）之内 ⇒ **不跨内含子边界，纯编码区事件**。
⚠ 这比 1126 文档中凭印象写的 `195,542,946-195,542,954` **左端差 1 bp**，本文的坐标以 CIGAR 为准。

**读框后果**【测】：8 不是 3 的倍数 ⇒ 移码。该区间覆盖**密码子 24 的末位 + 密码子 25、26 整个 + 密码子 27 的首位**。
重建该等位 CDS 并翻译：

| | 野生型 | Sample_5 第二条等位 |
| :--- | :---: | :---: |
| CDS 长度 | 615 bp | 607 bp |
| 蛋白长度 | **205 aa** | **45 aa**（终止于 PTC） |
| 前 24 残基 | `MAASTASHRPIKGILKNKTSTTSS` | **相同** |
| 第 25 位起 | `MVASAEQ...` | `VGRTAPRECRRGAEQKIPEVG*` |

⚠ 缺失虽然吃掉了密码子 24 的第 3 位，**残基 24 仍是 Ser** —— 重组后的密码子因 Ser 的 6 重简并性
恰好仍编码 Ser【测】。所以"前 24 残基相同"与"缺失覆盖密码子 24 末位"并不矛盾，
**移码的蛋白层面起点是残基 25**。

⇒ **结论（注释层面）**：Sample_5 两条等位分别为
**① Δ202（起始密码子 + 密码子 1–26 缺失）** 与 **② 8 bp 移码（45 aa 截短）**，
**两条都经过切割与 NHEJ 修复，且两条都预测为功能丧失**。
**凝胶看不出这个区别**【景】：8 bp 的长度差在为 202 bp 设计的胶上通常不分开。

**三条边界，缺一就读过头**：
1. PTC 落在 **exon 1**，距下游多个 exon-exon junction 很远 ⇒ **NMD 可能发生**，但本轮无 RNA 数据，**未评估**；
2. 以上全部是**注释层面的预测**，不是蛋白实测（无 Western、无功能试验）；
3. 前 24 残基保留意味着 `KGILK` 仍在该截短肽上 ⇒ **"短肽完全无害"也是未经验证的假设**。

### 5.4 缺失序列原文（供逐碱基比对，全部【测】）

> 全部由 `samtools faidx` 从 GRCh38 直读，复算命令见 §8。
> **已交叉核对**：`GATK.GRCh38/Homo_sapiens_assembly38.fasta` 与
> `human_gencode_v45/GRCh38.primary_assembly.genome.fa` 在本区间 `cmp` 逐字节相同。
>
> ⚠ **读之前先认链**：***PPP1R2*** **在负链**。IGV / 参考 FASTA 给的是**正链**；
> **mRNA sense = 负链 = 坐标从大到小读**。下面每段都给两条链，混用会得出完全错误的密码子。

#### 5.4.1 Δ202：`chr3:195,542,948–195,543,149`（202 bp）

**正链（参考 / IGV 方向），坐标升序**
```
195,542,948  CACCATAGAGGAAGTCGTAGAGGTCTTGTTCTTCAAGATCCCCTTGATGGGCCGGTGCGA
195,543,008  GGCCGTCGAGGCCGCCATTGCCGGGCGCTCCGGCTGTCGGCTCAGGGTCGCTGCTTGGCG
195,543,068  TGGGGTCCGCGAAGAGAAGGGTCGGCACAGCAGAGACTCGCAGGCAGCCGCAGATCCCGC
195,543,128  TCAGGGCTAAAGCGGCCGCAAC
```

**负链 = mRNA sense，5′→3′，坐标降序**
```
195,543,149  GTTGCGGCCGCTTTAGCCCTGAGCGGGATCTGCGGCTGCCTGCGAGTCTCTGCTGTGCCG
195,543,089  ACCCTTCTCTTCGCGGACCCCACGCCAAGCAGCGACCCTGAGCCGACAGCCGGAGCGCCC
195,543,029  GGCAATGGCGGCCTCGACGGCCTCGCACCGGCCCATCAAGGGGATCTTGAAGAACAAGAC
195,542,969  CTCTACGACTTCCTCTATGGTG
```

**同一 202 nt 按功能切成三段（mRNA sense）** —— 脚本断言三段首尾相接 == 上面整块，返回 `True`：

| 段 | 坐标 | 长度 | 序列（mRNA sense，5′→3′） |
| :--- | :---: | :---: | :---: |
| **[A] 5′UTR 部分** | `195,543,026–195,543,149` | 124 nt | `GTTGCGGCCGCTTTAGCCCTGAGCGGGATCTGCGGCTGCCTGCGAGTCTCTGCTGTGCCG` + `ACCCTTCTCTTCGCGGACCCCACGCCAAGCAGCGACCCTGAGCCGACAGCCGGAGCGCCC` + `GGCA` |
| **[B] 起始密码子** | `195,543,023–195,543,025` | 3 nt | `ATG` |
| **[C] 密码子 2–26** | `195,542,948–195,543,022` | 75 nt | `GCGGCCTCGACGGCCTCGCACCGGCCCATCAAGGGGATCTTGAAGAACAAGACCTCTACG` + `ACTTCCTCTATGGTG` |

124 + 3 + 75 = 202 ✓；**[B] + [C] = 78 nt = 恰好 26 个完整密码子** ⇒ 残余下游 CDS 仍在原读框上。
翻译 [B]+[C] = **`MAASTASHRPIKGILKNKTSTTSSMV`（26 aa）**，与 §5.1 所列一致。

⭐ **这一段序列顺带提供了一个独立互证**：翻译出的第 12–16 残基正是 **`KGILK`**，
与 UniProt **P41236**（人 I-2）的 N 端逐字相符 ⇒ **GENCODE 坐标、链向、读框三者同时被验证**，
而不是只靠注释文件自洽。它也直接说明 **Δ202 删掉的就是 §4.2 那个 PP1 对接基序本身**。

**逐密码子（可直接对着 IGV 点位核）**
```
codon  1 ATG M  195,543,023-195,543,025     codon 14 ATC I  195,542,984-195,542,986
codon  2 GCG A  195,543,020-195,543,022     codon 15 TTG L  195,542,981-195,542,983
codon  3 GCC A  195,543,017-195,543,019     codon 16 AAG K  195,542,978-195,542,980
codon  4 TCG S  195,543,014-195,543,016     codon 17 AAC N  195,542,975-195,542,977
codon  5 ACG T  195,543,011-195,543,013     codon 18 AAG K  195,542,972-195,542,974
codon  6 GCC A  195,543,008-195,543,010     codon 19 ACC T  195,542,969-195,542,971
codon  7 TCG S  195,543,005-195,543,007     codon 20 TCT S  195,542,966-195,542,968
codon  8 CAC H  195,543,002-195,543,004     codon 21 ACG T  195,542,963-195,542,965
codon  9 CGG R  195,542,999-195,543,001     codon 22 ACT T  195,542,960-195,542,962
codon 10 CCC P  195,542,996-195,542,998     codon 23 TCC S  195,542,957-195,542,959
codon 11 ATC I  195,542,993-195,542,995     codon 24 TCT S  195,542,954-195,542,956
codon 12 AAG K  195,542,990-195,542,992     codon 25 ATG M  195,542,951-195,542,953
codon 13 GGG G  195,542,987-195,542,989     codon 26 GTG V  195,542,948-195,542,950
```

**Δ202 的接合处（mRNA sense，两侧各 40 nt）** —— 左侧 40 nt 两型相同，右侧自第 2 位起分叉：
```
缺失型:  TGTCGTGACAACCGCTGCGGTAGCCGTTTCCGAGGCAGCA | GCGTCGGCCGAACAGCCCCGCGGGAATGTCGACGAGGAGC
野生型:  TGTCGTGACAACCGCTGCGGTAGCCGTTTCCGAGGCAGCA | GTTGCGGCCGCTTTAGCCCTGAGCGGGATCTGCGGCTGCC
```

#### 5.4.2 Sample_5 第二条等位：`chr3:195,542,947–195,542,954`（8 bp）

| 链 | 序列 |
| :--- | :---: |
| 正链（参考 / IGV） | `CCACCATA` |
| 负链（mRNA sense） | `TATGGTGG` |

**逐碱基的密码子归属**（负链；密码子 *n* 占 `195,543,026-3n` … `195,543,028-3n`）：

| 坐标 | 碱基（mRNA sense） | 归属 |
| :--- | :---: | :---: |
| `195,542,954` | `T` | 密码子 24 第 3 位 |
| `195,542,953` | `A` | 密码子 25 第 1 位 |
| `195,542,952` | `T` | 密码子 25 第 2 位 |
| `195,542,951` | `G` | 密码子 25 第 3 位 |
| `195,542,950` | `G` | 密码子 26 第 1 位 |
| `195,542,949` | `T` | 密码子 26 第 2 位 |
| `195,542,948` | `G` | 密码子 26 第 3 位 |
| `195,542,947` | `G` | 密码子 27 第 1 位 |

⇒ 密码子 24 末位 + 密码子 25、26 整个 + 密码子 27 首位 = 8 nt，**非 3 的倍数 ⇒ 移码**（§5.3）。

**接合处（mRNA sense，50 | 8 | 50）**
```
野生型:  CACCGGCCCATCAAGGGGATCTTGAAGAACAAGACCTCTACGACTTCCTC [TATGGTGG] CGTCGGCCGAACAGCCCCGCGGGAATGTCGACGAGGAGCTGAGGTGAGCG
缺失型:  CACCGGCCCATCAAGGGGATCTTGAAGAACAAGACCTCTACGACTTCCTC           CGTCGGCCGAACAGCCCCGCGGGAATGTCGACGAGGAGCTGAGGTGAGCG
```

⚠ **这段右侧 50 nt 已跨出 exon 1 的 CDS 末端（`195,542,904`），含内含子序列。**
**不要拿它直接数下游密码子** —— 读框后果必须在拼接后的 mRNA 上算（§5.3 的 45 aa 即如此得出）。

**§5.3 那个"看起来矛盾"的点在序列级闭合**：缺失吃掉密码子 24 的第 3 位（`T` @ `195,542,954`），
重组后该密码子 = `TC` + 下游首个保留碱基 `C` = **`TCC`**，仍编码 Ser（原为 `TCT`，亦 Ser）。
⇒ "前 24 残基与野生型相同"与"缺失覆盖密码子 24 末位"**不矛盾**，移码的蛋白层面起点是**残基 25**。

---

## 6. 这个基因的生物学对 P46 分析口径的影响

### 6.1 它不是读出指标，所以**不要**在报告里把它写成研究对象

客户的目的是定点制造双链断裂【客】。
⇒ 报告描述靶点时，主体是"**那 202 bp 的结构**"，**不是"PPP1R2 功能被破坏了"**。
后者即便成立，也不是客户问的问题，而且（§5.2）本轮数据支持不到那个层级。

### 6.2 ⚠ 但它可能是成像体系的混杂因素（需要在 Limitations 里留一句）

§4.3 若成立（I-2 参与分裂期调控、可能激活 Aurora A），那么**编辑克隆与对照克隆在 PPP1R2 剂量上本就不同**：
7 个编辑克隆中有 3 个（Sample_5/6/7 的先验）为双等位改动。
⇒ 任何"编辑组 vs 对照组"的表型差异，**除了已知的核型不同源（chr12p）之外，还多一条可能来源**。
**这不改变现有结论，只是 Limitations 需要点到。**

### 6.3 克隆存活性的一个未查项

单等位改动的克隆能长出来不奇怪；**双等位改动也能长出来**说明 PPP1R2 在 RPE-1 中**不是必需基因**【推】。
⇒ 这与 DepMap 类数据的预期是否一致，**本文未查**，列为 §9 待办。
**现在能说的只有**："本轮获得了双等位改动的存活克隆"这个**观察事实**。

---

## 7. 生物信息学：这个位点的四个实操特征

> 这一节是本文对分析最实用的部分 —— **这个位点有一个家族，而家族会咬人。**

### 7.1 家族构成：2 个旁系同源 + 6 个加工假基因（【测】，GENCODE v45）

| 基因 | 染色体 | 坐标 | 类型 |
| :--- | :---: | :---: | :---: |
| **PPP1R2** | chr3 | 195,514,428–195,543,386 | protein_coding（**本靶点**） |
| PPP1R2B | chr5 | 156,850,295–156,852,528 | protein_coding |
| PPP1R2C | chrX | 42,777,366–42,778,163 | protein_coding |
| PPP1R2P1 | chr6 | 32,879,216–32,879,833 | transcribed_processed_pseudogene |
| PPP1R2P2 | chr21 | 35,887,195–35,887,807 | processed_pseudogene（落在 RUNX1 区间内） |
| PPP1R2P4 | chr13 | 46,437,394–46,437,988 | processed_pseudogene |
| PPP1R2P5 | chr2 | 106,940,880–106,941,482 | processed_pseudogene |
| PPP1R2P6 | chr7 | 140,292,752–140,293,339 | processed_pseudogene |
| PPP1R2P10 | chr13 | 63,811,014–63,811,608 | processed_pseudogene |

**6 个假基因都是 processed pseudogene，长度 588–618 bp** ⇒ 它们是 mRNA 的反转录拷贝，
**所以同源落在 exon 序列上，而两个切点就在 exon 1 上。**

### 7.2 同源有多强：实测 76.8%–95.4%（【测】，blastn）

以**缺失区 202 bp ± 100 bp 侧翼**为 query：

| 位点 | identity | 比对长度 | bitscore |
| :--- | :---: | :---: | :---: |
| **PPP1R2B** (chr5) | **95.4%** | 346 bp | 553 |
| PPP1R2P1 (chr6) | 91.2% | 306 bp | 427 |
| PPP1R2P4 (chr13) | 87.1% | 348 bp | 408 |
| PPP1R2P10 (chr13) | 85.4% | 350 bp | 385 |
| PPP1R2P6 (chr7) | 85.1% | 316 bp | 349 |
| PPP1R2P2 (chr21) | 84.7% | 347 bp | 377 |
| PPP1R2P5 (chr2) | 76.8% | 354 bp | 279 |
| PPP1R2C (chrX) | 77.9% | 95 bp | 76 |

⇒ **这就是 P46 Stage 2c 的成对切点枚举"全基因组只有 6 处、且 5 处全是 PPP1R2 家族"的原因**【测·见 research record】。
**不是巧合，是旁系同源的必然结果** —— 两条 gRNA 的近似位点在家族成员上**同时**出现、**间距还与靶点相仿**，
因为那里就是同一段序列的拷贝。

### 7.3 ⚠ 但高同源 **没有** 变成比对歧义 —— 这一条与直觉相反，必须实测才知道（【测】）

`Ctr_1.markdup.bam`，排除 secondary/supplementary（`-F 0x900`）后的 MAPQ 分布：

| 窗口 | reads | MAPQ<20 | MAPQ=0 |
| :--- | :---: | :---: | :---: |
| **靶点 202 bp ± 100 bp** | 177 | **1.1%** | 0.6% |
| PPP1R2B（chr5，95.4% 同源） | 1,009 | 1.0% | 1.0% |
| PPP1R2P1（chr6） | 547 | 0.4% | 0.4% |
| PPP1R2P4（chr13） | 536 | 2.1% | 1.9% |
| 中性对照 chr7:50.0–50.1 Mb | 35,648 | 0.4% | 0.3% |
| 中性对照 chr2:100.0–100.1 Mb | 38,751 | 8.0% | 7.5% |

**判读**：85–95% 的同源**在读段层面不造成定位困难** —— 150 bp 的双端读段上，
5%–15% 的差异意味着每条读段带 **7–22 个判别性错配**，足够唯一定位。
靶点窗口（1.1%）与干净的中性区（0.4%）同量级，**远好于** chr2 那个随手选的"对照"窗口（8.0%）。

⚠ **两条方法学教训**：
1. **"有同源 ⇒ 必有比对歧义"是错的推理**，取决于同源程度与读长的关系。本例实测否定了它；
2. **"中性对照窗口"不能随手选** —— 我第一次选的 chr3:195.4–195.5 Mb（靶点附近 100 kb）
   MAPQ<20 高达 **26.9%**，因为 3q 末端本身是重复富集区。
   **拿它作对照会得出"靶点区异常干净"的假结论。** 对照窗口必须先自证干净。

⇒ **对报告的含义**：靶点判读不需要为同源性打折扣；
但 Stage 2b **分开记录 primary / secondary / supplementary** 的做法仍然必要 ——
它守的是"**定位歧义的证据本身**"，而非"歧义一定存在"。

### 7.4 高 GC / CpG island 区：对建库与覆盖的影响（【测】）

| 区段 | 长度 | GC | CpG o/e |
| :--- | :---: | :---: | :---: |
| **删除的 202 bp** | 202 bp | **65.8%** | 0.92 |
| 保留的 5′UTR + TSS | 176 bp | 69.3% | 1.24 |
| exon 1 全长 | 422 bp | 67.8% | 1.09 |
| TSS 上游 500 bp | 500 bp | 60.4% | 0.94 |
| 保留的 exon 1 CDS | 44 bp | 70.5% | 1.32 |

⇒ **整个靶点区符合 CpG island 判据**（Takai & Jones 阈值：GC ≥ 55% 且 CpG o/e ≥ 0.65）【测】+【推】。
⚠ 这是按序列组成算出来的判定，**未核对 UCSC 的 CpG island 注释轨**；要写进对外文字须先对一次注释。
**两条实际影响**：
1. **高 GC 区在 PCR 建库中系统性覆盖偏低** ⇒ 靶点窗口的绝对深度**不能直接与全基因组均深比**，
   深度类判据（`d_norm`）必须**用同一样本的同窗口做分母**，这正是 Stage 1 的做法【测·见 research record】；
2. **实测覆盖并不差**：靶点窗口 143–177 reads（Ctr_1 / Sample_7）⇒ 高 GC 没有造成判读盲区。

---

## 8. 可复算命令（全部为本文数字的来源）

```bash
G=/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/gencode.v45.annotation.gtf
FA=/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa
B=/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121/bam

# § 3  基因与转录本结构
awk -F'\t' '$9 ~ /gene_name "PPP1R2"/ && $3=="gene"' $G
awk -F'\t' '$9 ~ /ENST00000618156/ && ($3=="exon"||$3=="CDS"||$3=="start_codon")' $G

# § 7.1  家族成员
awk -F'\t' '$3=="gene" && $9 ~ /gene_name "PPP1R2(B|C|P[0-9]+)?"/' $G

# § 7.3  靶点区比对质量（与中性对照对比）
samtools view -F 0x900 $B/Ctr_1.markdup.bam chr3:195542848-195543249 \
  | awk '{t++; if($5<20) lo++} END{printf "reads=%d MAPQ<20=%.1f%%\n",t,100*lo/t}'
samtools view -F 0x900 $B/Ctr_1.markdup.bam chr7:50000000-50100000 \
  | awk '{t++; if($5<20) lo++} END{printf "reads=%d MAPQ<20=%.1f%%\n",t,100*lo/t}'

# § 5.3  Sample_5 第二条等位的 8 bp 缺失坐标（CIGAR 直读）
samtools view $B/Sample_5.markdup.bam chr3:195542900-195543000 \
  | awk '$6 ~ /8D/ {print $4,$6}' | sort | uniq -c | sort -rn | head

# § 5.4  两段缺失的参考序列（正链；负链序列 = 其反向互补）
samtools faidx $FA chr3:195542948-195543149   # Delta-202，202 bp
samtools faidx $FA chr3:195542947-195542954   # Sample_5 第二条等位，8 bp
samtools faidx $FA chr3:195543026-195543149   # 其中的 5'UTR 段，124 bp
samtools faidx $FA chr3:195543023-195543025   # 起始密码子（负链读作 ATG）

# § 5.4  两个 FASTA 给出同一序列的交叉核对（本文已跑，结果 identical）
FA2=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta
diff <(samtools faidx $FA  chr3:195542948-195543149 | tail -n +2) \
     <(samtools faidx $FA2 chr3:195542948-195543149 | tail -n +2) && echo identical
```

§5.2 / §5.3 / §5.4 的读框、翻译与序列切分计算用的一次性脚本留在 scratchpad
（`/tmp/claude-1001/-home-gao/76ac2071-2f9a-4347-a1de-4d68c31a7935/scratchpad/`，
文件名 `ppp1r2_anno.py`、`frame.py`、`s5_allele.py`、`del202_mrna.py`、`cpg.py`、`delseq.py`）。
**scratchpad 会被清理** ⇒ 若这些结论要进交付物，需把脚本移进 `scripts/` 并编号。

---

## 9. 待办 / 待核清单

| # | 要做什么 | 为什么要紧 | 判据 |
| :--- | :---: | :---: | :---: |
| 1 | **查 `KGILK` 基序与 Thr72 的原始文献**（§4.2） | 这是本文唯一进入"蛋白功能"层面的论断，目前凭印象 | 找到结构/生化原文则引用并统一残基编号；找不到则在对外文字中删去基序名称，只写"N 端 26 残基" |
| 2 | **查 I-2 的分裂期功能与 Aurora A 关系**（§4.3） | 决定 §6.2 那句 Limitations 怎么写 | 有文献 ⇒ 写成"已知参与分裂期调控"；无 ⇒ 只写"该基因功能未在本轮评估" |
| 3 | **查 PPP1R2 在 DepMap 中的必需性**（§6.3） | 支撑"双等位克隆能存活"这个观察的解释 | 非必需 ⇒ 与观察一致；若报为必需 ⇒ 须解释为何克隆能长出来 |
| 4 | **把 §5.3 的结论回填进 `dual_cut_per_cut_genotyping_1126.md` §5** | 那份文档明确写着"本文未做这个核对" ⇒ 不回填就会被后续会话再问一次 | 该节的 ⚠ 块改为指向本文 §5.3，并修正 8 bp 缺失坐标的左端（946 → 947） |
| 5 | **决定 §5.2 / §5.3 是否进客户报告** | 这是**新增的注释层面推断**，报告 1127 版已在内部审阅中 | 需用户拍板：是作为补充写进**报告的** §7.1 / §8.3，还是只留内部（注意本文自己也有 §7.1，此处指交付报告） |
| 6 | **决定报告 header 的 Tissue / Cell 怎么写**（§2A.4） | 现写 `hTERT RPE-1`（第 ① 层），只读 header 的人不会知道克隆带着整合的报告器；送样单 `cell_line` 字段原文是 `RPE-1_53BP1-FFR-mCherry monoclonal` | 需用户拍板：沿用现写法、改为送样单原文、或两者并列。**报告在内部审阅中，不擅自改** |

---

## 10. 与其他文档的分工

| 文档 | 管什么 |
| :--- | :---: |
| **本文** | **① PPP1R2 这个基因与蛋白**：结构、家族、靶点在基因结构中的位置、两条等位的注释层面后果、该位点的生物信息学特征；**② 材料与样本的命名**（§2A）：三层谱系、token 含义、误读点 |
| `sLK_of_CRISPR_Editing/dual_cut_per_cut_genotyping_1126.md` | **逐切点基因型判定的方法学**：为什么不能合并两切点求和、相位与克隆纯度怎么分开判 |
| `sLK_of_CRISPR_Editing/hTERT_RPE1_cell_line_biology_1125.md` | **细胞系的生物学**：核型（chr10q / chr12p）、p53 状态及其对分析口径的约束 |
| `sLK_of_CRISPR_Editing/mCherry_BP1_2_reporter_biology_1125.md` | **报告器构建体**：53BP1-FFR 的原理、元件坐标、同源屏蔽的代价 |
| `sLK_of_CRISPR_Editing/P31_P46_caged_dual_cut_wgs_rationale_1125.md` | **实验设计**：caged / dual-cut / WGS 为什么这么搭；对照非 CRISPR-naive 的论证 |
| `../scripts/0_sample_table.csv` | **样本名与送样单元数据的 SSOT** —— 本文 §2A 不复制这张表 |
| `../docs/P46_research_record_1125.md` | P46 各 stage 的做法、实测结果与产物路径 |

**⚠ §2A 与另两份文档的边界（避免同一信息存两处）**：
**本文 §2A 只管"名字指什么"**；
**细胞系的生物学性质**（核型、p53、必需性）归 `hTERT_RPE1_cell_line_biology_1125.md`；
**报告器构建体的内部结构**（元件坐标、屏蔽区）归 `mCherry_BP1_2_reporter_biology_1125.md`。
§2A 引用这两份的结论时只给指针，不重述论证。

**不一致时**：基因/转录本坐标与读框计算以**本文 §3 / §5** 为准（本文是首次按 GENCODE v45 逐条核的）；
Sample_5 的 8 bp 缺失坐标以**本文 §5.3** 为准（CIGAR 直读，修正了 1126 文档的左端）；
**样本名与送样单元数据以 `scripts/0_sample_table.csv` 为准**（本文 §2A 只解释名字，不存数据）；
各 stage 的分析结果以 research record 为准。
