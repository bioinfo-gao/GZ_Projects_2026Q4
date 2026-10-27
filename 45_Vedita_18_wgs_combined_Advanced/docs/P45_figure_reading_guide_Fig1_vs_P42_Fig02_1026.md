# 读图指南：P45 Fig1（locus overview）、P42 Fig02（DLD1 chr11 depth）与 P45 Fig5 接头碱基视图

> 内部文档，2026-10-26。对象：需要解读或向客户解释这两张图的我方人员。
> 本文只解释"图怎么读"，不修改任何交付产物；文末列出已发现的呈现问题与待决定的修改方案。
> 每条数值主张后给出来源文件与可复跑的核对命令（§10）。

## 0. 两张图在哪、看什么

| 图 | 文件 | 生成脚本 | 回答的问题 |
| :--- | :---: | :---: | :---: |
| P42 Fig02 | `/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis/custom_research_report_20261024/figures/P42_Fig02_DLD1_chr11_depth.pdf` | `42_.../scripts/20261023/10_make_delivery_figure_suite.py:113-129` | 有没有深度下降、范围多大（121.5 kb 宽视野） |
| P45 Fig1 | `/home/gao/projects_2026Q4/45_Vedita_18_wgs_combined_Advanced/custom_research_report_20261025/figures/Fig1_locus_overview.pdf` | `45_.../scripts/20261025/11_figures.py:78-140`；深度来自 `scripts/20261025/from_P44/9c_depth_boundary.py` | 相对对照降了多少、断点在哪、断点附近有无接头（junction）证据（37 kb 窄视野） |

两图用**同一批 BAM**（`projects_2026M7/41_Vedict_18_wgs/analysis/bam/S0x.markdup.bam`），看的是同一现象：
**S04/S05（+Dox 克隆）自 chr11 约 114.605 Mb 起测序深度约为对照的一半，并一直延续到 P42 窗口右端 114.675 Mb。**
两图外观不同几乎全部来自"纵轴除以什么"（§1），不是数据矛盾。

推荐读法：先看 P42 Fig02 确认"有、范围多大"，再看 P45 Fig1 读"降多少、断在哪、有无接头证据"。

---

## 1. 核心差异：两图的纵轴不是同一个量

| 维度 | P42 Fig02 | P45 Fig1 中间 panel |
| :--- | :---: | :---: |
| 纵轴 | 样本深度 ÷ **该样本全基因组平均深度** | 样本深度 ÷ 该样本**本地侧翼中位数**（chr11:114.40–114.58 Mb），再 ÷ **mean(S01–S03)** |
| 含义 | 每个样本的绝对相对深度，各自保留原貌 | +Dox 相对 no-Dox 对照的比值；对照共有的起伏被除掉 |
| bin 大小 | 1 kb | 200 bp |
| x 范围 | 114.554–114.675 Mb | 114.595–114.632 Mb |
| samtools depth 过滤 | `-Q 20 -G 2048 -s`（排除 supplementary，抑制 mate 重叠双计） | `-Q 20`，默认 flag（含 supplementary，不抑制重叠） |
| 画出的样本 | S01–S05 | S02、S03、S04、S05（S01 只进分母，见 §5） |

两点直接后果：

- **P42 能看到对照自身的形状**（repeat 区鼓包、上游组间差），P45 看不到——它们被分母抵消了。
- **P45 的比值在对照深度≈0 处无定义**，代码将其设为 NaN 不画（§3）。

过滤差异（`-G 2048 -s`）只影响深度绝对值的小幅偏移，不改变任何区段的相对关系。

### 1.1 用区段均值把两图对上

P42 数据（MAPQ（mapping quality，比对质量）≥ 20，1 kb bin），各区段均值：

| 样本 | 上游 114.554–114.600 | 114.602 凹陷 | repeat 内 114.606–114.625 | repeat 下游 114.635–114.675 |
| :--- | :---: | :---: | :---: | :---: |
| S01 pop | 1.19 (17.6×) | 0.05 | 1.32 (19.6×) | 1.08 (16.0×) |
| S02 ctrl | 1.18 (15.5×) | 0.02 | 1.49 (19.6×) | 1.10 (14.5×) |
| S03 ctrl | 1.20 (17.5×) | 0.03 | 1.61 (23.5×) | 1.10 (16.1×) |
| S04 +Dox | 1.42 (21.1×) | 0.04 | 0.77 (11.4×) | 0.59 (8.7×) |
| S05 +Dox | 1.54 (26.4×) | 0.01 | 0.81 (13.9×) | 0.53 (9.2×) |

括号内为原始深度。+Dox ÷ 对照均值：repeat 内 ≈ 0.52，repeat 下游 ≈ 0.54——与 P45 阶梯模型右侧比值 S04 0.431、S05 0.366 同量级（P45 的分母与过滤不同，数值不必逐位相等）。

---

## 2. 读 P42 Fig02 时必须知道的三个现象

### 2.1 上游区 +Dox 虚线高于对照（约 1.45 vs 1.2）

- 来源：S04/S05 在这一片的**原始深度本来就高**（21–26× vs 15.5–17.6×）。
- P45 用同染色体本地侧翼归一化后，该偏移消失（左侧比值 S04 0.985、S05 0.826）。
- **这是真实拷贝数差异还是文库覆盖偏好，未检查，两个方向都不下结论。**

### 2.2 对照在 114.603–114.626 鼓起到 1.3–1.6 ⇒ +Dox 看似"降了两次"

- 鼓包区正是 42 bp × 522 拷贝的串联重复 array（§4）。对照在此本身深度偏高。
- +Dox / 对照比值在 repeat 内（≈0.52）与 repeat 下游（≈0.54）基本相同。
- ⇒ P42 中 114.627 附近的"第二次下降"是 **repeat 区自身的深度膨胀造成的视觉错觉，不是第二个断点**。P45 除以对照后是一条平台。

### 2.3 两处所有样本都接近 0 的尖凹陷

见 §3。

---

## 3. 两处深度≈0：测序盲区，不是拷贝数变化

### 3.1 实测特征

| | 凹陷 1 | 凹陷 2 |
| :--- | :---: | :---: |
| 零深度区段（S02，MAPQ ≥ 0，逐碱基） | 114,601,509–114,602,941（中间有几处数碱基至数十碱基的非零小岛） | 114,632,287–114,632,935 |
| 参考序列 GC | 8%（任一 100 bp 窗口最高 17%） | 7%（最高 22%） |
| 参考序列样貌 | `TATATATTATATATTATAATTATATATAAT…` | `…AATATATAATATATAAATATTATATATTAC…` |
| 参考中 N | 无 | 无 |
| MAPQ ≥ 0 平均深度，18 样本 / 4 细胞系 | 0.00–0.21× | 0.00–0.33× |

### 3.2 判定

- **不是 mappability（可比对性）问题**：MAPQ ≥ 0（不过滤）也同样为 0 ⇒ 没有 read 落在这里，而不是 read 落下后被判为多重比对。
- **不是缺失**（推断）：DLD1、GM12878、Mutu I、DG75 四个无关细胞系同时为 0；共享缺失极不可能，而极端 AT + (TA)n 回文序列在建库 PCR / cluster 扩增中系统性丢失是 Illumina 的已知偏好。**两侧边界有无 soft-clip 断点未检查**，故"不是缺失"为推断。
- **对主结论的影响**：S04/S05 的深度阶梯（114.605 Mb）与两凹陷均不重叠。但凹陷 1 右缘（114,602,941）距 S04 断点阴影左缘（114,605,000）仅约 2.1 kb，断点定位不能向左延伸进这片盲区。

### 3.3 为什么 P42 是 V 字、P45 是空窗、凹陷 2 在 P45 里看不到

| | P42 Fig02 | P45 Fig1 |
| :--- | :---: | :---: |
| 凹陷 1 | 1 kb bin 114,602,000–114,602,999 几乎全零 ⇒ 值 0.01–0.05；折线连接相邻 bin ⇒ 画成 V | 6 个 200 bp bin（114,601,600–114,602,799）三对照均 ≤ 0.04；`11_figures.py:105` 将分母 ≤ 0.05 的 bin 设为 NaN ⇒ 断线，形成约 1.2 kb 空窗 |
| 凹陷 2 | 落在 114,632,000 bin，同样画成 V | 首个 bin 中心 114,632,100 > x 轴上限 `X1 = 114,632,000` ⇒ 在图外 |
| 分母会否为 0 | 不会（全基因组均值 14–17×） | 会（对照局部深度≈0） |

- P45 的空窗处理在逻辑上正确：对照无数据时"相对对照的比值"无定义；画成 0 反而会让人以为 +Dox 有缺失。
- 空窗两端 114,601,400 与 114,602,800 两个 bin 的对照均值仅约 0.10–0.14，刚过 0.05 门槛被画出，**比值极不稳定，不应解读**。
- P45 报告图注已写明 "the gap at 114.6015–114.603 Mb has no MAPQ ≥ 20 reads in any sample"，但未说明原因；P42 图未作任何标注。

---

## 4. P45 Fig1 顶部 panel：蓝条与橙条

### 4.1 画法

- 蓝条：分析目标区（target）chr11:114,604,000–114,625,500（`p45lib.py:19`）。
- 橙条：UCSC simpleRepeat（TRF，Tandem Repeats Finder 检出的串联重复）中 **period（重复单元长度）≥ 20 bp** 的全部记录，每条一个矩形（`11_figures.py:84-91`）。重叠/首尾相接的矩形连成一片，未被覆盖处成为缺口。
- 下方文字"tandem repeat array, 42 bp unit × 522 copies (114,603,768–114,625,660)"**只描述最长那一条记录**（`max(..., key=长度)`），不是整条橙色。

### 4.2 合并后的实际区间

| 区间（GRCh38，1-based） | 长度 | 内容 | 图上表现 |
| :--- | :---: | :---: | :---: |
| 114,601,509–114,602,359 | 851 bp | AT 富集串联重复（23/27 bp 单元，如 `TATATAATATATATTATAGATTA`） | 橙色第 1 段 |
| 114,602,360–114,602,433 | 74 bp | 只有 period 2/7/18 的 (TA)n 类记录，被 ≥ 20 阈值滤掉 | **前部小缺口** |
| 114,602,434–114,603,053 | 620 bp | AT 富集串联重复（18–163 bp 多种单元） | 橙色第 2 段 |
| 114,603,054 | 1 bp | 两条记录首尾相差 1 bp | 若见细缝，是相邻无边框矩形的抗锯齿接缝 |
| 114,603,055–114,625,660 | 22,606 bp | 47 bp 变体（至 114,604,208）+ 42 bp 主 array（521.9 拷贝），内部叠有大量 21 bp 半单元记录 | 橙色主体 |
| 114,625,661–114,631,681 | 6,021 bp | TRF **无任何记录**（任何 period） | **后部大缺口** |
| 114,631,682–114,632,022 | 341 bp | **同一 42 bp 家族**的残余拷贝（7.5 拷贝，86% 一致） | 右端小橙块，被 x 上限截断 |

### 4.3 读图要点

- **橙色起点 114.6015 ≠ 42 bp array 起点 114.6038。** 前两段是另一类重复：GC 约 8%，与 42 bp 单元（`GGTAACCACTGTTACCTGGTGGATAATAAGTATTGCCTCGTG`，GC 19/42 ≈ 45%）无可识别同源。凹陷 1（§3）几乎完全落在这两段之内（零深度区段止于 114,602,941，橙色第 2 段止于 114,603,053，相差约 110 bp）。
- **橙色本身没有画错**——它们确是 TRF 串联重复；**误导来自呈现**：单一颜色、无图例、标注文字置于整条橙色下方、报告图注用单数 "the tandem repeat array"。
- **后部 6 kb 缺口是判断"深度下降是否只是 repeat 比对假象"的关键区段**：在 TRF 未注释为串联重复的这 6 kb 里，S04/S05 比值仍约 0.4–0.5（P42 中延续到 114.675 Mb 仍约 0.55×）。这是观察，不等于证明了缺失；机制仍需接头证据。该段在 RepeatMasker 中是否含分散重复（LINE/SINE 等）未检查。
- **右端残余拷贝与主 array 同家族**：来自主 array 的 read 可能比对到这里（反之亦然），解读此附近的接头时需考虑。

---

## 5. P45 Fig1 中间 panel：深度比值

### 5.1 线与参考

- 紫 S04、蓝 S05：+Dox 克隆。左侧 ≈ 1，右侧 ≈ 0.43 / 0.37；log2 step ≈ −1.19 / −1.18，相当于约一半拷贝。
- 棕色细线 S02、S03：no-Dox 克隆。
- 两条点线：1.0 与 0.5 参考。

### 5.2 阴影 = 断点 95% CI（confidence interval，置信区间）

两段式阶梯模型拟合 log2 比值（114.590–114.640 Mb），1,000 次 block bootstrap（block = 5 bins，seed 44）。阴影范围为 `ci95_lo` 至 `ci95_hi + 200`（含最后一个 bin 的宽度）：

| 样本 | 断点 bin 起点 | 阴影范围 | 左侧比值 | 右侧比值 | log2 step |
| :--- | :---: | :---: | :---: | :---: | :---: |
| S04 | 114,605,200 | 114,605,000–114,605,600 | 0.985 | 0.431 | −1.193 |
| S05 | 114,605,800 | 114,605,600–114,606,200 | 0.826 | 0.366 | −1.176 |

两阴影首尾相接、不重叠 ⇒ 两克隆的断点位置在 200 bp 分辨率下可区分。两者都位于 42 bp 主 array 起点之后约 1.2–2.4 kb，即**在 repeat 内部**。

### 5.3 为什么没有 S01；为什么棕线不是独立对照

- S01 **参与了分母**（`11_figures.py:102`），但只画 S02/S03（`11_figures.py:103-111`）。不画的原因项目文档无记录；推断：
  1. **数学冗余**：分母是三者均值 ⇒ 每个 bin 中 S01/base + S02/base + S03/base ≡ 3，S01 线 = 3 − S02 线 − S03 线，画出无新信息。
  2. 图例按"克隆对克隆"组织（S04/S05 vs S02/S03），S01 是亲本群体。
- 同一原因使**棕线在 1 附近波动部分是构造决定的**，不能当作"对照平坦"的独立证据。对照有无阶梯请读阶梯表：

| 样本 | 左侧比值 | 右侧比值 | log2 step | CI 宽度 |
| :--- | :---: | :---: | :---: | :---: |
| S01 pop | 0.987 | 0.887 | −0.154 | 约 38 kb（无法定位） |
| S02 ctrl | 0.982 | 1.022 | +0.058 | 约 48 kb |
| S03 ctrl | 0.975 | 1.069 | +0.133 | 约 44 kb |

  这只说明三者之间无一明显偏离；真正独立的阴性比较是 GM12878 组（S07–S10 对 S06，log2 step −0.28 至 −0.69；S09/S10 的 CI 较窄，位于 114.627–114.638 Mb 即 repeat 下游一带，S07/S08 的 CI 宽约 46 kb 无法定位；该组另需单独解读）。
- S01 原始深度可在 P42 Fig02 直接看（深蓝实线，§1.1 表）。

---

## 6. P45 Fig1 底部 panel：接头事件

**这个 panel 回答：S04/S05 深度下降的地方，有没有 read 直接跨过一个接头（junction）？**
答案：有，但证据有限。258 个点里绝大多数是 repeat 内比对歧义造成的 1-fragment 背景；真正要读的只有 114.605–114.607 Mb 的 3 个三角形——两个 +Dox 克隆各有接头，且与各自的深度断点位置吻合。

### 6.1 编码

| 元素 | 含义 | 来源 |
| :--- | :---: | :---: |
| 每个点 | `junction_table_target.tsv` 中一个目标区接头候选 | 共 527 条；0 fragment 的 269 条（Not_supported）不画 ⇒ 258 个点，全部在 x 范围内 |
| x | 接头落在目标窗口一侧的断端坐标 | `w_side()`，`11_figures.py:65-68` |
| y | read rescue 后 **S01–S10 任一样本**中支持该接头的最大独立 DNA fragment 数 | 列 `max_after_frag_S01_S10` |
| 颜色 | 证据分级（§6.2） | 列 `classification` |
| ▽ / △ | 仅 S04 / 仅 S05 有，S01–S03 均无 | 列 `DLD1_dox_comparison` 以 `case_specific` 开头（含 underpowered），共 39 个 |
| ○ | 其余全部 | 含 GM12878 独有、DLD1 对照独有、多细胞系共有 |

y 的计数口径（报告 §5）：一端 read 须在接头两侧各跨 ≥ 20 bp，且对接头模板的比对得分（AS，alignment score）比其全基因组最佳比对高 ≥ 10，才算一次支持；fragment = 独立 DNA 分子，同一对 read 的两端只计一次。

### 6.2 分级标准（报告 §5）

| 分级 | 标准 |
| :--- | :---: |
| Well supported | 单样本 ≥ 3 fragment；≥ 2 种证据（split read、discordant pair、组装 contig）；两侧 100 bp 侧翼均可唯一放置；断点不确定度 ≤ 10 bp |
| Limited | ≥ 2 fragment，侧翼唯一，不确定度 ≤ 10 bp；但 < 3 fragment 或仅一种证据 |
| Ambiguous repeat | ≥ 2 fragment，但某侧翼或伙伴端可同样好地放到多个位置，或断点可在相同序列内移动 > 10 bp |
| Weak | 1 fragment；列出，不解读 |

### 6.3 读图顺序

**第一步：把 y = 1 的那一排当背景。**
254 个 Weak 点从 114.604 铺到 114.626，宽度与 42 bp array 几乎一致。在 22 kb 串联重复里，单个 read 很容易被放到错误的拷贝上，表现为"接头"。它们不是 254 个结构事件。

**第二步：记住点不全是 DLD1 的。**
258 个点中 145 个只出现在 GM12878，只在 DLD1 出现的为 100 个；y 也取 S01–S10 最大值，一个圆点的高度可能来自 GM12878。**只有三角形可读为"DLD1 +Dox 特异"。**

**第三步：只看 y ≥ 2 的 4 个点。**

| 图上的点 | junction | target 侧断端 | 内容 | 分级 / fragment | 样本 | 解读 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| y=4 紫灰 ▽ | J0922 | 114,605,218 | chr11 → (TTAGGG)n 端粒重复 | Ambiguous_repeat / 4 | S04 | chr11 端唯一；伙伴端在 chrX、chr2、chr4 末端及 chr10 间隙端粒样序列上同样好放置，表中 chrX:156,030,577 只是代表性放置之一；标准计数 0 → rescue 后 4；**落在 S04 断点阴影内** |
| y=2 红 △ | J0359 | 114,605,973 ↔ 114,606,349 | 376 bp fold-back 倒位 | Limited / 2 | S05 | **落在 S05 断点阴影内** |
| y=3 红 △ | J0786 | 114,607,192 → chr21:20,648,241 | 染色体间接头 | Limited / 3 | S05 | 位于 S05 断点下游约 1.2 kb，**不在阴影内**；DELLY、GRIDSS、SvABA 三个 caller 均检出（即旧报告 EV0623） |
| y=2 红 ○ | J1036 | 114,609,273 → 114,609,283 | 10 bp 局部事件 | Limited / 2 | **S03（对照）** | 对照自身事件，与 +Dox 无关；**最易被误读为 +Dox 证据** |

前三个即报告 §7.4 的"三个 +Dox 接头"。

**第四步：与中间 panel 的阴影对齐。**
J0922 在 S04 阴影（114,605,000–114,605,600）内，J0359 在 S05 阴影（114,605,600–114,606,200）内 ⇒ 每个 +Dox 克隆各有一个接头恰在自己深度下降处，空间上自洽。这就是报告中 "coherent with the copy-number step" 的含义；该判定还要求第二个条件——断点下游原杂合位点在该样本中 ≥ 90% 变纯合（见 Fig2 BAF，B-allele frequency）。

**第五步：看哪里没有点。**
- 两克隆无任何共享接头（S04 的点在 S05 中不存在，反之亦然）⇒ 两次断裂是独立事件。
- 图例中的 "Well supported" 在本窗口内一个都没有（图例按全部分级生成）。

### 6.4 能说与不能说

| 可以说 | 不能说 |
| :--- | :---: |
| S04、S05 各有一个 2–4 fragment 的接头，位置与各自深度断点吻合 | "断点已确认"——没有独立实验；对同一数据重比对/重计数不是独立验证 |
| J0922 的 chr11 端唯一，伙伴端是端粒重复 | "J0922 连到 chrX"——该坐标只是代表性放置 |
| 254 个 Weak 点主要反映 repeat 内比对歧义 | "Weak 是假的"——单 fragment 不足以判真假，只是不解读 |
| J0786 是三 caller 支持的 S05 接头 | "J0786 解释了 S05 的深度下降"——它在 CI 之外 |

本数据能达到的最高结论层级是 **"limited, coherent with the copy-number step"**；直接确认需跨接头 PCR（报告 §8.4 已设计引物）。

---

## 7. Fig5 接头碱基视图（`figures/junction_views/`）

### 7.1 位置与范围

- 目录：`/home/gao/projects_2026Q4/45_Vedita_18_wgs_combined_Advanced/custom_research_report_20261025/figures/junction_views/`，`Fig5_<junction_id>.{png,pdf}`，**17 个接头 × 2 格式 = 34 个文件**。报告正文只嵌入 `Fig5_J0922`（报告图注 "Figure 5 (one of 17 …)"）。
- 入选规则（`11_figures.py:245`）：rescue 后 S01–S10 任一样本 ≥ 2 fragment 的接头 = 目标区表 4 个（§6.3 表）+ 窗口其余区段表（`junction_table_window_context.tsv`）13 个。

| junction | 位置 | 分级 / fragment | 出现样本 | 与 +Dox 的关系 |
| :--- | :---: | :---: | :---: | :---: |
| **J0922** | chr11:114,605,218 → 端粒重复 | Ambiguous_repeat / 4 | S04 | +Dox 特异，在 S04 断点 CI 内 |
| **J0359** | 114,605,973 ↔ 114,606,349 倒位 | Limited / 2 | S05 | +Dox 特异，在 S05 断点 CI 内 |
| **J0786** | 114,607,192 → chr21 | Limited / 3 | S05 | +Dox 特异，CI 外 |
| **J0591** | 114,602,528 → chr16 | Limited / 2 | S04 | +Dox 特异但 underpowered |
| **J0902** | 114,629,858 → chrX | Ambiguous_repeat / 2 | S05 | +Dox 特异，repeat 下游 |
| J1036 | 114,609,273 局部倒位 | Limited / 2 | S03 | 对照自身事件 |
| J1365、J0533 | 114,590,721 / 114,590,799 | Well_supported / 10、8 | S06–S10 | GM12878 全体共有，细胞系特征 |
| J1592 | 114,555,232 → chrX PAR1 | Limited / 2 | S07（另见 Mutu I、DG75） | 无关；解读见 §7.5 |
| 其余 8 个（J0771、J0907、J1635、J1775、J2143、J0715、J1219、J1474） | 窗口各处 | Limited / Ambiguous | 多样本共有或单个 GM12878 克隆 | 与 baseline 共有或 GM12878 克隆事件 |

### 7.2 图的结构（以 J0922 为例）

1. **标题行**：`J0922 BND A=chr11:114605218:+ B=chrX:156030577:- (representative placement shown) class=… microhomology=GT inserted=-`
   - `:+`/`:-` 为断端方向：`chr11:114605218:+` = 沿 chr11 正链读到该位置后跳入伙伴序列。
   - `microhomology`：两侧共享的碱基，可归任一侧 ⇒ 断点只能定到区间（J0922 为 114,605,216–218）。`inserted`：两端之间不属于任一侧的插入碱基。
2. **第 2 行**：颜色说明 + 伙伴端状态（"B placement unique"，或 "AMBIGUOUS PARTNER … has N alternative placements besides the coordinate shown (N+1 equivalent in total)"——2026-10-26 修正后的措辞，见 §9 第 7 条）。注意 "placement unique" 指**接头模板 100 bp 侧翼**在参考中唯一，不等于支持 read 本身可唯一定位（§7.5）。
3. **junction 参考行**：断点左右各 70 bp 的拼接序列。酒红 = A 侧（chr11），深青 = 微同源/插入，蓝 = B 侧（伙伴）。竖虚线 = 接头，位于微同源右边界。
4. **read 行**：标签 "样本 Fxx/Rn"。同一 F 号的 R1/R2 是同一分子，**只计 1 个 fragment**（J0922：7 行 = F01–F04 = 4 fragment）。`.` = 与接头序列相同；字母 = 错配（显示 read 实际碱基）；`-` = 缺失；空白 = 未覆盖。只画**通过判定**的 read 端（两侧各连续跨 ≥ 20 bp 且接头模板 AS 比全基因组最佳高 ≥ 10），故有的 fragment 只出现 R1 或 R2。所有 read 已转到接头模板方向。

### 7.3 四个判据

1. **每行连续穿过虚线**，不在虚线处截断 ⇒ split read 直接证据。
2. **虚线两侧的 `.` 越长越好**（≥ 20 bp 是门槛；J0922 右侧约 40–80 bp）。
3. **错配看分布**：
   - 同一分子 R1/R2 在**不同**位置的错配 = 随机测序错误（J0922 F04：R1 为 C/C/T，R2 为 G/G/G，位置互不重合；AS 接头 139/135 vs 全基因组 83/74）。
   - **多个独立分子在相同位置有相同错配** = read 的真实来源序列与参考不同，相关一端的定位可疑（J1592，§7.5）。
4. **独立分子起点不同** ⇒ 不是 PCR duplicate（去重已排除）。

可以说："N 个独立分子在碱基层面直接连接 A 与 B"。不能说："已验证"（同一数据的重比对不是独立实验）；伙伴端标为 representative 时，不能说"连到该坐标"。

### 7.4 为什么有的 R1/R2 两行完全相同、有的 fragment 只有单端（以 J0922 为例）

**结论**：
- **F02、F03 的 R1/R2 完全相同**：这两个 DNA 分子只有 87 bp、94 bp，短于 150 bp 读长。R1 和 R2 都从分子一端读到另一端，再读进 Illumina 接头（adapter）。adapter 在比对时被 soft-clip，图上不画 ⇒ 两条 read 覆盖同一段插入序列，画成两行一样的点。
- **F01 只有 R2**：F01 的 R1 止于断点前 1–3 bp，未跨过接头，不满足"两侧各 ≥ 20 bp"，因此不画。R2 跨过接头进入端粒 23 bp。

**实测**（`key_junction_read_evidence.tsv` 的 primary 记录；序列转回测序方向后搜索 Illumina adapter `AGATCGGAAGAGC`）：

| fragment | 读端 | 比对（primary） | adapter 在 read 中的起点 | 插入片段 | 是否跨接头 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| F01 | R1 | chr11:114,605,066 `150M` | 无 | 约 176 bp | **否**：止于 114,605,215，接头在 114,605,216–218 |
| F01 | R2 | chr11:114,605,092 `127M23S` | 无 | 同上 | 是；端粒侧 23 bp，刚过 20 bp 门槛 |
| F02 | R1 | chr11:114,605,153 `63S66M21S` | 第 87 位 | **87 bp** | 是 |
| F02 | R2 | chr11:114,605,153 `66M84S` | 第 87 位 | **87 bp** | 是 |
| F03 | R1 | chr11:114,605,153 `56S66M28S` | 第 94 位 | **94 bp** | 是 |
| F03 | R2 | chr11:114,605,153 `66M84S` | 第 94 位 | **94 bp** | 是 |
| F04 | R1 | chr11:114,605,131 `88M62S` | 无 | 约 165 bp | 是 |
| F04 | R2 | chrX 端粒 `71S79M`（supplementary：chr11:114,605,146 `73M77S`） | 无 | 同上 | 是 |

**从 CIGAR 读出短插入片段**（F02 = 66 bp chr11 + 21 bp 端粒 = 87 bp）：
- R2（正链，flag 163）`66M84S`：66 bp 比对 chr11；84 bp soft-clip = 21 bp 端粒（TTAGGG，不属于 chr11）+ 63 bp adapter。
- R1（反链，flag 83）`63S66M21S`（参考方向）：63 bp adapter（反向互补后在前）+ 同一段 66 bp chr11 + 同一段 21 bp 端粒。
- 两者起点同为 114,605,153，chr11 部分同为 66 bp，跨接头后同为 21 bp ⇒ 转到接头模板方向后完全重合。F03 同理（66 + 28 = 94 bp）。
- **F04 对比**：插入约 165 bp > 读长；R1 起于 114,605,131、R2 起于 114,605,146，重叠约 135 bp，两端各错开约 15 bp ⇒ 图上两行长度不同。§7.3 判据 3 "R1/R2 错配位置不同 = 测序错误" 正是在这段重叠区里成立：同一碱基被两条 read 各读一次而结果不一致。

**对计数的意义**：
- **R1/R2 相同 ≠ 两份证据**：同一分子，fragment 只计 1（F02 = 1，F03 = 1）。J0922 的 4 fragment = F01–F04 四个分子；图上 7 行 read 并未被当作 7 份证据。
- 短插入片段两端读到同一段序列 = 同一碱基测了两遍：提高碱基准确度，不增加分子数。
- **F01 的 R1 是锚点**：证明该分子的 chr11 一侧确在此处，但本身未跨接头，故不画。F01 的支持完全依赖 R2 进入端粒的 23 bp——4 个 fragment 中跨接头最短的一个。
- F02、F03 只有 87、94 bp，明显短于文库典型插入长度。adapter 已剪除、参与比对的只是真实插入序列，故不影响判定；但它意味着 J0922 一半的支持来自很短的 DNA 片段，每个在端粒侧只有 21–28 bp。

### 7.5 案例：`Fig5_J1592`——看起来 Limited，实际不应按克隆事件解读

**结论**：J1592 是 GM12878 克隆 S07 的接头候选，在目标区左侧约 49 kb（P45 Fig1 x 范围之外），与 DLD1/+Dox 无关。图上两个信号要求降低信任：(1) 支持 read 来自三个互不相关的细胞系；(2) 所有 read 在相同位置带相同错配。分级 Limited 只因 S07 恰有 2 fragment 且模板侧翼可唯一放置。报告正文未提及 J1592。

**标题**：`A=chr11:114555232:-  B=chrX:1919999:-  microhomology=-  inserted=AGTATAAATAATGCACATAATGCAG`
- B 端在 chrX PAR1（拟常染色体区）。
- 无微同源，**25 bp 插入**（深青段）。

**参考行**：酒红侧末尾 `CTCAAAAAAAAAA`——A 侧位于 **AluSc8（SINE）**，这串 polyA 是 Alu 3′ 端尾巴；蓝侧位于 **L1MA7（LINE）**；A 侧另有 chr17 上 94% 相似的片段重复（segmental duplication）。结构为"Alu 尾 → 25 bp → L1"，两端均在可移动元件内，是比对歧义与群体多态的高发结构。

**read 行**：

| 行 | 细胞系 | 说明 |
| :--- | :---: | :---: |
| S07 F01/R1 | GM12878 克隆 5 | 仅 R1 通过 |
| S07 F02/R1、R2 | GM12878 克隆 5 | 同一分子，计 1 fragment |
| S11 F03/R1、R2 | **Mutu I 群体** | 背景细胞系 |
| S15 F04/R1、R2 | **DG75 群体** | 左侧未覆盖或 soft-clip |
| S17 F05/R1 | **DG75 克隆 13** | 背景细胞系 |

**关键：共享错配。** 酒红侧第 9 位 `C`、第 15–16 位 `GG`，蓝侧约 45 bp 处 `T`，更远处 `AC`——**4 个细胞系的 5 个独立分子在完全相同位置有完全相同的错配**（对比 J0922 的随机错配）。推断：这些 read 来自一份与 chr11 这段 Alu 高度相似但不完全相同的拷贝（另一 Alu 位点，或参考未收录的多态 Alu–L1 结构）；chr11 只是参考中最接近的位置。佐证：原始比对 MAPQ 均为 0，XA 列出几十个跨染色体替代位置；接头模板 AS 92–125 比全基因组最佳 55–74 高 20 以上而通过规则——**该规则只证明 read 最像接头模板，不保证接头模板就是其真实来源**。

**样本层面**：

| 样本 | 支持 fragment | 参考 fragment | 状态 |
| :--- | :---: | :---: | :---: |
| S01–S05（DLD1） | 0 | 8–21 | absent（S01、S03 参考 < 10，underpowered） |
| S06（GM12878 群体，S07 亲本） | 0 | 7 | **underpowered** |
| S07（GM12878 克隆 5） | 2 | 16 | present，AF 0.11 |
| S08–S10 | 0 | 10–11 | absent |
| S11–S18（背景） | 3 个样本有支持 | — | Mutu I、DG75 |

- 表中标签 "S07: case_specific_underpowered"：相对 S06 像是 S07 独有，但 S06 只有 7 个参考 fragment，不能排除 S06 也有。
- 克隆性杂合事件预期 AF ≈ 0.5，实测 0.11。
- 四个 SV caller 均未检出（`n_callers=0`）。
- 叠加跨三个无关细胞系出现 ⇒ "S07 克隆特异事件"不成立；更符合"带共同序列的多态/参考未收录结构"或"重复序列系统性比对伪像"，两者本图无法区分。

| 可以说 | 不能说 |
| :--- | :---: |
| S07、S11、S15、S17 中有 read 连续跨过一个"Alu polyA → 25 bp → L1"接头 | S07 有克隆特异的 chr11–chrX 易位 |
| 与 DLD1/+Dox 无关 | chr11 端确在 114,555,232（read 与参考有共享差异，来源存疑） |
| 支持 read 带一致的序列差异，原始 MAPQ 为 0 | "这是假阳性"（也可能是真实多态结构，本数据不下结论） |

---

## 8. 一句话总结

先用 P42 回答"有无、范围"，再用 P45 回答"降多少、断在哪、附近有无接头"。两图一致：约 50% 深度下降，起点 114.605–114.606 Mb，位于 42 bp repeat array 内部；两处深度≈0 的凹陷是 AT 富集测序盲区，与处理无关；断点附近只有 1–2 个 +Dox 特异、2–4 fragment 的接头，证据有限且受 repeat 干扰。

---

## 9. 已发现的呈现问题（第 7 条已修正；其余未修改，待用户决定）

| # | 图 | 问题 | 可选修改 |
| :--- | :---: | :---: | :---: |
| 1 | P42 Fig02 | repeat 区对照鼓包导致"两次下降"错觉；无对照比值 | 图注说明，或加一个 +Dox / mean(对照) 比值 panel |
| 2 | P42 Fig02 | 两处 V 形凹陷未标注，可能被读成缺失 | 灰色阴影 + "AT-rich, no coverage in all 18 samples" |
| 3 | P45 Fig1 顶部 | AT 富集段与 42 bp array 同色相连；标注只对应最长记录；图注单数 | 方案 A（推荐）：AT 富集段改灰/斜线并标注，42 bp 家族保留橙色；B：只画 42 bp 家族；C：仅改图注 |
| 4 | P45 Fig1 中间 | 未说明为何不画 S01；棕线非独立对照 | 图注补 "S01 enters the baseline but is not drawn" |
| 5 | P45 Fig1 中间 | 空窗原因未说明；凹陷 2 在 x 范围外 | 图注补原因；或 x 上限外扩至 114,633,500 |
| 6 | P45 Fig1 底部 | 圆点混有 GM12878 事件；图例含本窗口不存在的 Well supported | 图注说明，或按细胞系区分形状/透明度 |
| 7 | Fig5（伙伴端有歧义的 3 张：J0922、J1635、J1474） | 第 2 行 "has 3 equivalent placements; the coordinate shown is one of them" 与数据不符：代表位置 1 + 替代 3 = 4 处（J0922：chrX、chr10、chr2、chr4，报告正文亦为 4 处） | **✅ 已修正 2026-10-26**：改为 "has N alternative placement(s) besides the coordinate shown (N+1 equivalent in total)"；3 张图的 PNG/PDF 已替换（stage11 与交付目录），manifest 与 md5 已更新；旧版归档于 `superseded_for_comparison/analysis/v6_20261026_fig5_partner_placement_count_wording/` |
| 8 | Fig5_J1592 等 | 分级只看 S01–S10 最大值，未计入背景细胞系（S11–S18）中的出现；J1592 因此与 J0359 同为 Limited，正文无解释 | 图上标注 background presence；或把 ≥ 2 个无关细胞系出现的接头单独归类 |
| 9 | Fig5 全部 | read 标签只有样本号，看不出 S11/S15/S17 属于 Mutu I/DG75 | 标签加细胞系，或按细胞系着色 |

两图均已在交付目录中。任何修改前，须先将现有 PNG/PDF 与脚本快照归档到各自项目的 `superseded_for_comparison/`。

---

## 10. 核对命令

```bash
P42=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
P45=/home/gao/projects_2026Q4/45_Vedita_18_wgs_combined_Advanced
REF=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta
SAM=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/samtools

# §1 两边的 depth 命令与归一化
sed -n 18,21p $P42/scripts/20261023/3b_run_discovery_depth.sh
sed -n 1,10p  $P45/scripts/20261025/from_P44/9c_depth_boundary.py
sed -n 101,111p $P45/scripts/20261025/11_figures.py

# §1.1 区段均值的输入表
head -3 $P42/analysis/20261023_fullcohort_repeat_resolve/03_coverage/discovery_1kb_normalized.tsv

# §3 凹陷：18 样本 MAPQ≥0 深度
for d in $P42/analysis/20261023_fullcohort_repeat_resolve/03_coverage/discovery_per_base/*/; do
  zcat $d/inclusive_mapq0.depth.tsv.gz | awk -v n=$(basename $d) \
    '$2>=114601500&&$2<=114602900{a+=$3;c++} $2>=114632300&&$2<=114632900{b+=$3;e++} END{printf "%s %.2f %.2f\n",n,a/c,b/e}'; done
# §3 参考序列（检查 N 与 GC）
$SAM faidx $REF chr11:114601509-114602941 chr11:114632287-114632935
# §3.3 P45 空窗 bin
awk -F'\t' '$1==20 && $2>=114601000 && $2<=114603400' \
  $P45/analysis/20261025_inherited/P44_r1_stage9_supporting/depth_boundary/depth_bins_200bp.tsv | cut -f1-8

# §4 橙条记录（period≥20 为第 6 列）
zcat $P45/refs/simpleRepeat.txt.gz | awk -F'\t' '$2=="chr11" && $4>114595000 && $3<114632000 && $6>=20' | cut -f3-11,17

# §5 阶梯与 CI
grep -E "^20" $P45/analysis/20261025_inherited/P44_r1_stage9_supporting/depth_boundary/boundary_estimates.tsv

# §6 接头表（≥2 fragment）
awk -F'\t' 'NR==1||$26>=2' $P45/analysis/20261025_p45/stage10_tables/junction_table_target.tsv | cut -f1,3,4,6,26-28

# §7 Fig5：文件清单、17 个入选接头、单个接头全部字段与支持 read
ls $P45/custom_research_report_20261025/figures/junction_views/ | wc -l
for t in junction_table_target.tsv junction_table_window_context.tsv; do
  awk -F'\t' 'NR>1&&$26>=2' $P45/analysis/20261025_p45/stage10_tables/$t | cut -f1,3,4,6,26-28; done
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[i]=$i} $1=="J1592"{for(i=1;i<=NF;i++) if($i!="") print h[i]"="$i}' \
  $P45/analysis/20261025_p45/stage10_tables/junction_table_window_context.tsv
awk -F'\t' '$1=="J1592"{print $3,$4,$7,$13,$14,$18,$19}' \
  $P45/analysis/20261025_p45/stage10_tables/key_junction_read_evidence.tsv | sort -u
# §7.4 J0922 各 fragment 的 primary 比对与 adapter 位置（插入片段长度）
awk -F'\t' '$1=="J0922" && $10=="primary" && $7!="0"' $P45/analysis/20261025_p45/stage10_tables/key_junction_read_evidence.tsv \
  | awk -F'\t' '!s[$4 FS $8]++{print $4,"R"$8,$9,$11":"$12,$13,$NF}' \
  | python3 -c "
import sys
rc=lambda s:s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
for l in sys.stdin:
    fr,r,fl,pos,cig,s=l.split(); ro=rc(s) if int(fl)&16 else s
    print(fr,r,fl,pos,cig,'adapter_at',ro.find('AGATCGGAAGAGC'))"
```
