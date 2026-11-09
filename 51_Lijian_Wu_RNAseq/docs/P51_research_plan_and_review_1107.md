# P51 分析 — Lijian Wu 小鼠细胞系 bulk RNA-seq：control vs Crtc2 line × 4 处理

- **物种 / 材料**：Mus musculus；客户未提供名称的已建立小鼠细胞系（GRCm39 / GENCODE M35）
- **设计**：2 个细胞系（control、Crtc2 line）× 4 种处理（un、TNF、PGE2、adenosine）× 3 个生物学重复 = 24 样本；样本表 SSOT = [`scripts/0_sample_table.csv`](../scripts/0_sample_table.csv)
- **当前结果目录**：[`Data_Analysis_20261108/`](../Data_Analysis_20261108/)；交付目录 [`custom_research_report_20261108/`](../custom_research_report_20261108/)（**内部审阅中、未发送**）
- 旧版文档（含过程记录）已归档：[`superseded_for_comparison/docs/v1_20261108_process_log_with_data_issues/`](../superseded_for_comparison/docs/v1_20261108_process_log_with_data_issues/)

---

## 1. 核心结论

1. **Crtc2 line 的 Crtc2 mRNA 在 4 种条件下都下降，但没有消失**：降到 control 的 0.53–0.67 倍（log2FC −0.58 ~ −0.93）。修饰性质（KO / KD / 其它）客户未说明，所以只称 "Crtc2 line"。
2. **Crtc2 line 的主要表型是 CREB/cAMP 应答减弱**：经典 CREB/CRTC 靶基因（Nr4a2、Crem、Pde4b、Pde4d、Sik1）在 4 种条件下都显著更低。未处理时细胞系之间的差异与 PGE2 应答方向相反：同时属于两组 DEG 的 740 个基因里有 91% 方向相反。
3. **PGE2 应答的缺失最清楚**：control 对 PGE2 有强应答（2,174 DEG）；Crtc2 line 只保留约 58% 的应答幅度。PGE2 条件下两系差异的 2,105 个 DEG 中，66% 是 control 的 PGE2 应答基因，其中 97% 方向相反。也就是说，这一对比测到的主要是"Crtc2 line 没有做出的那部分 PGE2 应答"。
4. **TNF 条件下两系差异最大（5,905 DEG）**：Crtc2 line 仍对 TNF 有应答，但幅度较小（约 74%）。TNF 诱导的炎症基因明显更低，例如 Il6 −3.37、Ptgs2 −2.48、Cxcl10 −1.54。
5. **adenosine 在 control 中的作用很弱（254 DEG，|log2FC| ≥ 1 的只有 13 个），在 Crtc2 line 中几乎没有作用**：Crtc2-adenosine 与 Crtc2-un 的 PCA 质心距离为 0.34，与重复之间的距离相当。Crtc2-adenosine vs control-adenosine 的 DEG 主要来自两系本身的差异，其中 68% 也是未处理条件下的 DEG，且方向 100% 相同。
6. **通路层面，两系的差异在所有条件下方向一致**：Crtc2 line 的 TNF/NF-κB、inflammatory response、IFN-γ 程序更低，E2F/G2M 细胞周期程序更高（TNF 条件下细胞周期一项不显著）。这与 control 对 PGE2/adenosine 的应答方向正好相反。

---

## 2. 对比与方法

**对比**：只做客户送样单 3G 列出的 6 个（后者为 baseline）。

| # | 对比 | 回答的问题 |
| :--- | :---: | :---: |
| 1 | Crtc2-un vs control-un | 未处理时两系的差异 |
| 2 | Crtc2-TNF vs control-TNF | TNF 下两系的差异 |
| 3 | Crtc2-PGE2 vs control-PGE2 | PGE2 下两系的差异 |
| 4 | Crtc2-adenosine vs control-adenosine | adenosine 下两系的差异 |
| 5 | control-PGE2 vs control-un | control 的 PGE2 应答 |
| 6 | control-adenosine vs control-un | control 的 adenosine 应答 |

没做的对比：control-TNF vs control-un、Crtc2 各处理 vs Crtc2-un（客户未列）。下文涉及 TNF 应答的结论只基于 PCA 和组均值，属于描述性结果。

**方法**：

| 步骤 | 选择 |
| :--- | :---: |
| 上游 | nf-core/rnaseq 3.15.1，STAR（`--outFilterMultimapNmax 3`）+ Salmon gene-level estimated counts；strandedness 自动推断（24/24 为 reverse） |
| 基因类型过滤 | 去核糖体蛋白 / 线粒体 / 非编码与预测基因：57,132 → 26,584 |
| 低表达过滤 | 每个对比单独决定检验集合：本对比 6 个样本中 ≥3 个 ≥10 counts，或全体 ≥30% 样本 ≥10 counts；`independentFiltering=FALSE`，在集合内做 BH。全局并集 12,327 个基因，每个对比检验 11,823–12,172 个 |
| 模型 | 一个 DESeq2 `~ Group`（8 水平）模型，6 个对比都从中提取；Wald test；`lfcShrink(type="ashr")` |
| DEG 阈值 | padj ≤ 0.05 且 \|log2FC\| ≥ 0.263（1.2 倍）；同时报告 \|log2FC\| ≥ 1 的计数 |
| ORA | clusterProfiler GO BP/MF/CC + KEGG；UP / DOWN / ALL 分开做；背景 = 该对比的检验集合 |
| GSEA | 排序值 = DESeq2 未收缩的 Wald `stat`；GO BP、KEGG、MSigDB Hallmark（小鼠同源） |

阈值用 1.2 倍的理由：Crtc2 本身只变化约 −0.6 log2FC，大量方向一致的 CREB 靶基因变化也小于 2 倍，用 2 倍阈值会漏掉这一层主要信号。

---

## 3. 结果

### 3.1 PCA：PC1 = TNF 轴，PC2 = 细胞系轴

- PC1 解释 93% 的方差，PC2 解释 4%。组内平均距离 0.62，组间平均距离 20.58，比值 33.2。
- **PC1 由 TNF 驱动**：TNF 组的 PC1 均值为 27.9，其余组在 −12.1 到 −5.0 之间。
- **PC2 完全按细胞系分开**：所有 Crtc2 line 样本的 PC2 都在 0.1–6.4，所有 control 样本都在 −5.5 到 −1.7，两系没有重叠。
- **最近邻检查**：23/24 个样本的最近邻是同组重复。唯一的例外是 Crtc2-un-1 → Crtc2-adenosine-3，这两组本来就分不开（见下表），因此这不代表样本错标。

处理造成的位移（PCA 质心距离）：

| 处理 | control（处理 vs un） | Crtc2 line（处理 vs un） | Crtc2 / control |
| :--- | :---: | :---: | :---: |
| TNF | 47.36 | 32.96 | 0.70 |
| PGE2 | 10.55 | 4.64 | 0.44 |
| adenosine | 2.45 | 0.34 | 0.14 |

两系之间的距离：un 5.42，TNF 19.82，PGE2 10.63，adenosine 6.98。另外，Crtc2-PGE2 与 control-un 只相距 2.58，比 Crtc2-PGE2 与 Crtc2-un 的距离（4.64）还近。

### 3.2 DEG 数

| 对比 | 检验基因 | Up | Down | DEG 总数 | 占检验基因比例 | \|log2FC\| ≥ 1 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | 11,883 | 461 | 690 | **1,151** | 9.7% | 246（56 / 190） |
| Crtc2-TNF vs control-TNF | 12,172 | 2,804 | 3,101 | **5,905** | 48.5% | 1,353（556 / 797） |
| Crtc2-PGE2 vs control-PGE2 | 11,878 | 798 | 1,307 | **2,105** | 17.7% | 539（102 / 437） |
| Crtc2-adenosine vs control-adenosine | 11,823 | 495 | 809 | **1,304** | 11.0% | 305（53 / 252） |
| control-PGE2 vs control-un | 11,937 | 1,278 | 896 | **2,174** | 18.2% | 625（436 / 189） |
| control-adenosine vs control-un | 11,872 | 170 | 84 | **254** | 2.1% | 13（11 / 2） |

- 4 个 Crtc2 vs control 对比都以下调为主（下调占 52–62%）；放到 2 倍阈值后，下调占比升到 59–83%。
- **4 个 Crtc2 vs control 对比共享 656 个 DEG，其中 612 个方向一致**（404 个都下调，208 个都上调）。这是不随处理改变的细胞系核心差异。
- 未处理条件下差异最大的基因：下调方向有 Edil3 −2.90、Ramp3 −2.85、Aif1l −2.76、Wnt10b −2.65、Tmem119 −2.54、Rgs9 −2.50；上调方向有 Adrb2 +3.04、Lag3 +2.38、Plpp1 +1.87。

### 3.3 Crtc2 本身

DESeq2 normalized counts，3 个重复的均值：

| 条件 | control | Crtc2 line | Crtc2 / control | log2FC（对比 1–4） |
| :--- | :---: | :---: | :---: | :---: |
| un | 1,377 | 920 | 0.67 | −0.58* |
| TNF | 4,078 | 2,141 | 0.53 | −0.93* |
| PGE2 | 1,817 | 1,195 | 0.66 | −0.60* |
| adenosine | 1,433 | 943 | 0.66 | −0.59* |

- TNF 会上调 Crtc2 本身（组均值：control 3.0 倍，Crtc2 line 2.3 倍；TNF vs un 未做检验）。
- PGE2 在 control 中也上调 Crtc2（+0.38*）；adenosine 不改变 Crtc2（+0.01）。

### 3.4 标志基因面板

数值为 log2FC，`*` 表示显著，`n.t.` 表示表达太低、不在该对比的检验集合中。

| 基因 | 1 un | 2 TNF | 3 PGE2 | 4 adenosine | 5 control PGE2 应答 | 6 control adenosine 应答 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Crtc2 | −0.58* | −0.93* | −0.60* | −0.59* | +0.38* | +0.01 |
| Nr4a2 | −0.40* | −1.04* | −0.85* | −0.80* | +0.85* | −0.03 |
| Crem | −0.67* | −1.42* | −0.57* | −0.73* | +0.19 | −0.10 |
| Pde4b | −0.59* | −1.88* | −1.29* | −0.68* | +1.00* | −0.00 |
| Pde4d | −1.19* | −2.55* | −2.22* | −1.87* | +2.17* | +0.08 |
| Sik1 | −0.66* | −1.08* | −0.91* | −0.60* | +0.37* | −0.01 |
| Il6 | −0.24 | −3.37* | −2.37* | −0.54* | +2.69* | +0.01 |
| Ptgs2 | n.t. | −2.48* | −2.06* | n.t. | +2.33* | n.t. |
| Cxcl2 | n.t. | −0.63* | −0.98* | n.t. | +2.14* | n.t. |
| Cxcl10 | −0.70* | −1.54* | −1.46* | −0.19 | +0.58* | −0.01 |
| Ccl5 | −0.87* | −1.52* | −1.83* | −1.43* | +2.24* | +0.95* |
| Nfkbia | +0.10 | −1.17* | −0.39* | −0.05 | +0.67* | +0.27* |
| Tnfaip3 | +0.02 | −0.76* | −0.41* | −0.10 | +0.48* | +0.04 |
| Fos | +1.59* | +2.34* | +0.16 | +1.56* | +0.19 | +0.06 |
| Fosb | +1.35* | +3.00* | +0.39* | +1.75* | +0.05 | +0.01 |

如何读这张表：

- **CREB 靶基因**（Nr4a2 / Crem / Pde4b / Pde4d / Sik1）在 4 种条件下都显著更低。差距在 TNF 和 PGE2 下最大，与"PGE2 → cAMP → CREB/CRTC"这条通路受损一致。Pde4d 的变化最大且最稳定：PGE2 在 control 中把它上调 +2.17，而 Crtc2 line 在 PGE2 下比 control 低 −2.22。
- **炎症基因**（Il6 / Ptgs2 / Cxcl2）在 control 中受 PGE2 诱导，在 Crtc2 line 中 PGE2 和 TNF 条件下都明显更低。这些基因在未处理时差别不大，说明受损的是诱导应答，基础表达基本不变。
- **NF-κB 反馈基因**（Nfkbia / Tnfaip3）只在 TNF 下降低（−1.17 / −0.76），这与 Crtc2 line 的 TNF 应答幅度较小一致。
- **Fos / Fosb 方向相反**：在 un、TNF、adenosine 下都在 Crtc2 line 中更高，在 PGE2 下差别很小。这是面板里唯一与 CREB 靶基因方向相反的一组即刻早期基因（见 §4 推断 3）。

### 3.5 Crtc2 line 保留了多少处理应答（描述性）

方法：取 control 中有应答的基因，比较 Crtc2 line 的应答幅度（处理 vs un 的组均值 log2 比）与 control 的应答幅度，看两者之间的斜率。斜率 = 1 表示与 control 一样，0 表示完全没有应答。

| 处理 | 基因集 | 基因数 | 斜率（保留比例） | 相关系数 | 应答不足一半的基因占比 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| PGE2 | 对比 5 的 DEG | 2,174 | 58% | 0.86 | 33% |
| adenosine | 对比 6 的 DEG | 254 | 32% | 0.46 | 55% |
| TNF | control-TNF vs control-un 组均值差 ≥2 倍（未做检验） | 3,950 | 74% | 0.92 | 33% |

这张表与 PCA 位移比（TNF 0.70、PGE2 0.44、adenosine 0.14）给出相同的排序：TNF 保留得最多，adenosine 几乎全部丢失。

### 3.6 两系差异与处理应答的关系

| 交叉比较 | 两组都显著的基因 | 方向相反的比例 | 全基因 log2FC 相关 |
| :--- | :---: | :---: | :---: |
| 对比 1（un 下两系差异） × 对比 5（PGE2 应答） | 740 | 91% | −0.45 |
| 对比 1 × 对比 6（adenosine 应答） | 146 | 90% | −0.09 |
| 对比 3（PGE2 下两系差异） × 对比 5 | 1,384 | 97% | −0.77 |

- 对比 3 的 DEG 中 66% 属于对比 5 的 DEG ⇒ PGE2 下两系的差异，主要就是 Crtc2 line 没有做出的那部分 PGE2 应答。
- 对比 4 的 DEG 中 68% 也是对比 1 的 DEG，方向 100% 相同 ⇒ adenosine 下两系的差异基本等同于两系的基础差异，adenosine 本身几乎没有贡献。

### 3.7 通路

Hallmark GSEA 的 NES（正值表示在前一组更高；n.s. 表示不显著）：

| 对比 | TNFA_SIGNALING_VIA_NFKB | INFLAMMATORY_RESPONSE | INTERFERON_GAMMA_RESPONSE | E2F_TARGETS | G2M_CHECKPOINT |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | −1.56 | −2.15 | −1.65 | +2.00 | +1.86 |
| Crtc2-TNF vs control-TNF | −2.10 | −1.62 | −1.96 | n.s. | n.s. |
| Crtc2-PGE2 vs control-PGE2 | −2.21 | −2.22 | −2.12 | +2.10 | +1.97 |
| Crtc2-adenosine vs control-adenosine | −1.81 | −2.24 | −1.98 | +2.36 | +2.39 |
| control-PGE2 vs control-un | +2.59 | +2.19 | +2.35 | −3.31 | −3.15 |
| control-adenosine vs control-un | +2.75 | +2.17 | +2.46 | −3.79 | −3.53 |

ORA GO BP 每个方向的前几条：

| 对比 | 前一组更高 | 前一组更低 |
| :--- | :---: | :---: |
| Crtc2-un vs control-un | GPCR signaling; leukocyte / mononuclear cell migration | GPCR signaling; regulation of immune effector process; positive regulation of cell adhesion |
| Crtc2-TNF vs control-TNF | monocarboxylic acid / fatty acid metabolism; lipid biosynthesis | cytokine-mediated signaling; response to bacterium; regulation of immune effector process |
| Crtc2-PGE2 vs control-PGE2 | chromosome segregation; nuclear division | immune effector process; response to bacterium |
| Crtc2-adenosine vs control-adenosine | chromosome segregation; nuclear division | adaptive immune response; immune effector process |
| control-PGE2 vs control-un | immune effector process; cytokine-mediated signaling | chromosome segregation; nuclear division |
| control-adenosine vs control-un | T cell activation; cytokine-mediated signaling | chromosome segregation; sister chromatid segregation |

显著条目数（FDR < 0.05）：

| 对比 | GO BP up | GO BP down | KEGG up | KEGG down | GSEA GO BP | GSEA KEGG | GSEA Hallmark |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | 113 | 533 | 6 | 17 | 266 | 10 | 15 |
| Crtc2-TNF vs control-TNF | 137 | 516 | 14 | 18 | 75 | 24 | 14 |
| Crtc2-PGE2 vs control-PGE2 | 177 | 1,067 | 12 | 38 | 479 | 32 | 18 |
| Crtc2-adenosine vs control-adenosine | 186 | 604 | 26 | 28 | 465 | 4 | 18 |
| control-PGE2 vs control-un | 1,024 | 253 | 43 | 12 | 1,190 | 102 | 28 |
| control-adenosine vs control-un | 273 | 203 | 20 | 4 | 1,008 | 100 | 25 |

读法：
- control 的 PGE2 应答与 adenosine 应答在通路层面是同一个模式（炎症/IFN 程序上调，细胞周期下调），只是强度不同：基因层面 adenosine 的 DEG 少，GSEA 上却同样强。
- 这个模式在 Crtc2 line 中整体反转。TNF 下细胞周期一项不显著，原因可能是 TNF 本身对增殖的影响掩盖了两系差异（这是推断）。
- adenosine 应答 254 个 DEG，GSEA GO BP 却有 1,008 条：说明 adenosine 的作用是很多基因的小幅协同变化，单个基因很少过阈值。

---

## 4. 解读：实测事实与推断分开

**实测事实**（对应文件见 §6）：§1–§3 中的所有数字。

**我的推断**（下面均为假设，需要实验证实）：

1. **Crtc2 line 的表型符合 CRTC2 功能部分丧失**：Crtc2 mRNA 只降到 ~0.6 倍，但 CREB 靶基因与 PGE2（cAMP 上游）应答一起减弱，PGE2 应答受损（58%）远比 TNF（74%）严重。对一个只降低 mRNA 的细胞系来说，这样的表型强度提示：要么蛋白层面的损失比 mRNA 更彻底，要么是剪接或截短等改变造成的功能丧失。**RNA-seq 区分不了这两种情况**。
2. **adenosine 在 Crtc2 line 中几乎没有应答**：这与 adenosine 经 A2A/A2B → Gs → cAMP → CREB/CRTC 通路发挥作用相符。由于 control 本身的应答也很弱（254 DEG），"丢失 68%"这一数值的不确定性最大。
3. **Fos / Fosb 方向相反**：两者同样受 CREB 调控，却在 Crtc2 line 中更高。可能是 CRTC2 依赖性靶基因（Nr4a/Pde4/Sik1）与即刻早期基因的 CREB 依赖程度不同，或者是负反馈（Sik1、Pde4 降低 → cAMP/CREB 磷酸化代偿性升高）。这是推断，数据无法区分。
4. **TNF 应答减弱**：CRTC2 与 NF-κB 在炎症基因启动子上存在协同（Il6、Ptgs2 都是 CREB + NF-κB 双重调控的基因）。这可以解释 TNF 下差异最大的原因，但本数据不能证明直接作用。

**结论的边界**：
- §3.5 的保留比例是基于组均值的描述性统计，不是 genotype × treatment 交互检验；如需正式检验，可以在同一数据上增加 `~ line + treatment + line:treatment` 模型（客户未要求，属可选增量）。
- 通路名称只说明主题，不说明机制；mRNA 的变化也不等同于蛋白或功能的变化。
- n = 3：没进 DEG 列表不代表没有变化。

---

## 5. 建议的后续验证（写进报告 §8.4）

| 待回答的问题 | 实验 |
| :--- | :---: |
| CRTC2 蛋白是否缺失 | Western blot（control vs Crtc2 line） |
| 关键基因变化能否在独立实验中重现 | qPCR：Pde4d、Nr4a2、Il6、Crtc2（un / PGE2 / TNF） |
| Crtc2 line 对各处理的应答是否显著减弱 | 同一数据上做 genotype × treatment 交互检验（可选增量） |
| 通路差异是否有功能后果 | 增殖实验；IL-6 ELISA |

---

## 6. 文件索引（均在 `Data_Analysis_20261108/` 下）

| 内容 | 路径 |
| :--- | :---: |
| DEG 汇总 | [`DE_PCA_Results/DEG_summary_all_comparisons.csv`](../Data_Analysis_20261108/DE_PCA_Results/DEG_summary_all_comparisons.csv) |
| 每个对比的全表 / 显著表 / 火山图 / top50 热图 | `DE_PCA_Results/<comparison>/`（6 个） |
| PCA 与样本距离 | [`DE_PCA_Results/PCA_AllSamples_Overview.pdf`](../Data_Analysis_20261108/DE_PCA_Results/PCA_AllSamples_Overview.pdf)、`PCA_coordinates_overview.csv`、`SampleDistance_Heatmap.pdf` |
| 应答保留比例（§3.5） | [`DE_PCA_Results/Crtc2_line_treatment_response_retention.csv`](../Data_Analysis_20261108/DE_PCA_Results/Crtc2_line_treatment_response_retention.csv) |
| 标志基因面板（§3.4） | [`DE_PCA_Results/CREB_CRTC_target_and_inflammatory_marker_panel.csv`](../Data_Analysis_20261108/DE_PCA_Results/CREB_CRTC_target_and_inflammatory_marker_panel.csv) |
| 富集汇总 | [`Enrichment/Enrichment_summary_all_comparisons.csv`](../Data_Analysis_20261108/Enrichment/Enrichment_summary_all_comparisons.csv)；明细在 `Enrichment/<comparison>/{GO,KEGG,GSEA}/` |
| 客户报告（英文，审阅中） | [`Crtc2_Line_Treatment_RNAseq_Report_1108.md`](../Data_Analysis_20261108/Crtc2_Line_Treatment_RNAseq_Report_1108.md) |
| 脚本 | `scripts/4_run_DE_PCA.R`、`5_run_enrichment.R`、`6_generate_research_report.R`、`8_pre_delivery_audit.R`、`9_build_delivery.sh` |

§3.2 的 656 / 612 共享核心 DEG，以及 §3.6 的交叉比较数字，由 `*_ALL_genes.csv` 独立重算（2026-11-08），与报告数值一致。
