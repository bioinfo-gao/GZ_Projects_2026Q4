# P48 — LatchBio AI Biologist (Cancer Biology) take-home：研究计划与复盘

Plan Date（创建）: 2026-11-05
项目路径: `/home/gao/projects_2026Q4/48_Saul_LatchBio/`
状态: **计划阶段，未获授权运行正式分析**（用户未说"自动/全流程/全过程"）

## 更新记录

- 2026-11-05 — 创建。完成 Task 1 的数据勘察与陷阱验证；Task 2 仅有设计方向，数据源待定。
- 2026-11-05 — 据用户"BBB 通透性无统一标准"论点重写 §2.2：新增 2×2 列联表（**`(N,Y)`=0，
  72/72 满足逻辑包含**）、TableS1 内部三来源异质性、定量数值零命中的验证；
  **修正上一轮"零例反向⇒定义不同"的欠定推断**（§2.2.1）；新增 §2.3.1 grader 设计；
  §3.4 用同一论点解掉 Task 2 的数据源阻塞，两题合并为一条连贯主张。
- 2026-11-05 — 执行 §4 第 1–3 步：§2.1.4 判明 Bupropion 为**仓库真实录入错误**
  （标准答案改以附表为准，17 个冲突）；§1.1 确认**仓库无 license**；
  新增 §2.3.2（四篇判据，Doniger = CNS 活性口径**已证实**）与
  **§2.3.3（文献同时在用五个阈值，10%~200% 横跨 20 倍；阈值 0.3 恰等于实验误差；
  模型 RMSE 0.37 已超过该边距）**；§4.1 记录三项未闭合文献缺口。
- 2026-11-05 — **澄清一处我自己引入的错误假设：对方从未要求 JSON 格式。** `Email.txt:19` 只要求
  "挑一个熟悉的公开数据集、设计最多两道题、一个明确答案、一个只有专家才看出的陷阱、展示科学深度"；
  JSON schema 是我从其公开仓库反推的**内部格式**，对上它是加分信号但不是作答要求，
  而它把题目与 `ground_truth` 分置文件两头，**限时文本框里最难用**（用户指出）。
  ⇒ 渲染器新增 `--submit` 模式，按"题目 → 答案 → 陷阱 → 难点 → 测什么 → 答案来源 → 评分 →
  数据 → 局限"的线性顺序生成纯文本，答案紧跟题目；故意不含 `version_history`（内部迭代史）。
  产出 `Submission/submission_both.txt`（390 行）与两份单题版。**用户 2026-11-05 将其放入独立的 `Submission/`** —— 与 `docs/` 物理隔开，因为本文含大量不该给对方看的内容（v1–v4 为何被推翻等），考试当场不能开错文件。**由同一份 JSON 生成，不手写第三份副本。**
- 2026-11-05 — 两项可读性/诚实性修补：① 渲染器 `scripts/2_render_eval.py` 改为**按可见性排序**
  并加 `👁 模型可见` / `🔒 模型看不到` 两道分界 —— 起因是用户读中译时把 `notes.what_this_tests`
  误读成"题面都标注了，太简单"，而 `notes` 与 `ground_truth` 同在一个 JSON 里、本就是答案侧，
  模型只拿到 `task`；旧渲染把两侧排成一篇连续文档，必然反复制造这个误判。四份当前有效渲染已重排。
  ② Task 1 的 `notes.known_gaps` 补入**无法移除的部分捷径**：TableS3 文件名
  `...Based-on-CNS-Activity.xlsx` 把标注判据写在名字里，而题面必须给出该路径 ⇒ 仔细读文件名
  即可得到指向 `implication` 的提示。**结论是结果怎么读**：`relation` 答对不等于看出结构零，
  必须与五个整数计数一起读。`task` 与 `ground_truth` 的 sha1 实测未变（纯增量注释）。
- 2026-11-05 — **Task 2 再推翻重做为 v4**（新增 §3.7，§3.6 标记为已被取代）。用户追问
  「为什么要告诉 AI 这些事实」⇒ v3 题面那四条"临床事实"**全是 GBM 一般生物学而非患者特异信息**，
  其中第 1 条（增强区 BBB 已破坏）**恰恰就是让 BBB 成为诱骗项的那条知识，等于自己拆了自己的陷阱**。
  v4 整块删除事实清单与三药机制描述，抗 VEGF 选项缩为裸方向（原选项名里
  `unchanged_because_margin_barrier_already_intact` 本身也在陈述前提）。题面 2799 → 2186 字符。
  **这是同一错误的第三次（词形→解法→知识）**，根因是把"充分指定"误当成"前提喂到嘴里"；
  已把界线写成 v4 的 `notes.specification_principle`（§3.7.2）。v3 已逐字节归档。
- 2026-11-05 — **Task 2 推翻重做为 v3**（新增 §3.6，§3.3 标记为已被取代）。两条独立理由：
  ① 我自审发现 v2 六字段中四个白送/通识、决定性字段**纯排除法可解**，且题面那句"选三药共有属性"
  主动给排除路径立了路标；② 用户指出 **BBB 在脑瘤中本不重要（glioma 多伴 BBB 破坏），应作诱骗项**。
  由此得出可复用判据：**正确答案与文献主流叙事同向的陷阱不可能难**。v3 = 5 字段，BBB 作诱骗项
  + 两个方向反转字段（持续抗 VEGF 对暴露/对内在敏感性，与 vascular normalization 框架相反）
  + 一个自洽性检验字段；为可辩护性明确限定"超出正常化窗口之后"。v2 已逐字节归档
  （[v2_20261105_elimination_solvable/](../superseded_for_comparison/evals/v2_20261105_elimination_solvable/)，
  含中英两份渲染 + 缺陷说明）；新增 [evals/README.md](../evals/README.md) 声明哪个是当前有效版。
- 2026-11-05 — 新增 eval 渲染器 `scripts/2_render_eval.py` 与 `tmp/` 可读副本（中英各一套）；
  据用户选定方案 A 修掉 Task 2 v1 的**字面泄题**（`growth-arrested` 与选项 `quiescence` 近义、
  字段名含 `arrest` ⇒ 可纯词形匹配），v1 已归档。
- 2026-11-05 — 用户拍板 Task 1 用**蕴含关系框架**，并指示 §4.1 两项文献缺口暂时忽略。
  执行 §4 第 4、6 步：§3.2 两点 glioma 收紧经文献确认（另得"增强非暴露代理"、
  "bevacizumab 恢复屏障的悖论"、"低氧维持静止干细胞样群体"三条更好的陷阱）；
  §3.3 落定 Task 2 设计并写出两个 eval 文件 + 标准答案复现脚本（§4.2，共 7 个文件，
  均过校验）；§4.3 记录交付前必须做的双向验证与模型实测。
- 2026-11-05 — 新增 [scripts/2_render_eval.py](../scripts/2_render_eval.py)（eval 的人类友好渲染器，
  仅标准库；`--md` 产物自带"非真身"来源标记），按用户要求在 `tmp/` 生成 2 份可读副本。
  **渲染后发现 Task 2 存在字面泄题**：事实 #4 的 `growth-arrested` 与正确选项
  `target_cell_quiescence` 近义、字段名亦含 `arrest` ⇒ 可纯词形匹配答对。
  用户选定方案 A，已修为 v2（§3.5）；**v1 原件逐字节归档**于
  `superseded_for_comparison/evals/v1_20261105_lexical_leak_growth_arrested/`（含 README
  与可跑的并排 diff 命令）。机械复查：`arrest` 全题面 0 命中。

---

## §0 任务性质与硬约束

| 项 | 内容 | 来源 |
| :--- | :---: | :---: |
| 岗位 | AI Biologist - Cancer Biology (Applications), LatchBio | `Email.txt` |
| 第一轮形式 | HackerRank take-home，出题而非答题 | `Email.txt` |
| 时长 | **90 min，按 Continue 才计时，一次做完，不得刷新或后退** | `LatchBio.txt` |
| 邮件另称 | 2–3 小时 | `Email.txt` |
| 产出 | 最多 2 道题，须有唯一明确答案 + 只有专家能识破的陷阱 | `Email.txt` |
| 允许用 AI | Task creation 部分允许 | `Email.txt` |
| 第二部分 | coding / knowledge quiz / problem-solving，**不用 AI** | `Email.txt` |

**未解决的冲突**：90 min 与 2–3 h 不一致，且"能否在计时前准备好数据与答案"未确认。
已建议用户向 Saul 发邮件确认；**此问题的答复会改变本计划的可用范围**（若不允许预先准备，
下述验证产物不得带入考试）。

**无样本表**：本项目无测序样本，`scripts/0_sample_table.csv` 不适用。输入是公开数据仓库，
数据量记在 §1。

---

## §1 输入数据（实测，2026-11-05）

| 数据集 | 路径 | 体积 | 文件数 | 许可 |
| :--- | :---: | :---: | :---: | :---: |
| CASE-BBB-prediction-Data | `github.com/bioinformatics-gao/CASE-BBB-prediction-Data` | **11 MB** | **20** | **无 license（已确认，见 §1.1）** |

### 1.1 License 状态（2026-11-05 经 GitHub API 确认）

| 检查项 | 结果 |
| :--- | :---: |
| `license` 字段 | **None** |
| `/license` API 端点 | **HTTP 404** |
| 仓库树内 LICENSE / COPYING 文件 | **无**（唯一匹配项 `TableS2-High-Level-Group-**Terms**-...` 是假阳性） |
| 公开性 / forks / stars | public / 2 / 3 |
| 创建、最后 push | 2016-06-05、2016-10-16 |

**无 license = 默认版权保留**，第三方严格说并未获得再分发授权。
**但仓库 owner 账号 `bioinformatics-gao` 应为用户本人**（论文一作 + 账号名含 gao；
**待用户确认**）。⇒ 若确为用户所有，**用户自己加一个 license（CC0 / CC-BY / MIT）即可解除**，
这是一键可解的事，不是阻塞。提交题目前须完成。

- 勘察副本（可丢弃）：`/tmp/claude-1001/-home-gao--claude/3db4114b-4ef3-41e2-9441-cc56ac6ad44c/scratchpad/bbb_repo_dl/repo`
- 解析脚本：同一 scratchpad 下 `bbb_scripts/{scan2,overlap,inter,conflict}.py`
- **选它的理由**：公开可取 + 用户本人一作（*Bioinformatics* 2017, btw713），
  面试时可把来源与判定依据讲到任意深度。这是它优于 TCGA / DepMap 的地方。
- **⚠ license 未标注**：出题若要求受试者下载该仓库，需先确认可否被第三方 benchmark 引用。
  这是提交前必须解决的一项。

### 本项目不使用的数据

用户明确说明：**日常只分析病人与动物数据，公开数据挖掘仅用过 DrugBank。**
⇒ TCGA-PAAD / TCGA-GBM / Ivy GAP / DepMap / GSE115978 **全部排除**，
不满足"你熟悉的数据集"这一要求。此前我提过的黑色素瘤 scRNA-seq 方向已撤回
（经四个项目目录检索，用户无此项目）。

---

## §2 Task 1 — BBB 标签的证据层级（陷阱已验证）

### 2.1 实测事实

全部来自 `2) Intermediate-data-set/` 下真正喂给模型的两张特征矩阵（各 89 列，第 89 列为 BBB 标签）：

| 项 | Inter-Table4（训练+验证） | Inter-Table5（独立测试） |
| :--- | :---: | :---: |
| 药物数 | **213** | **161** |
| 标签分布 | Y 139 / N 74 | Y 76 / N 85 |
| 标签来源 | Li 2005 (138), Abraham 2006 (53), Subramanian & Kitchen 2003 (21), Wang 2015 (1) | Doniger 2002（附带 MW / volume / surface area） |

| 交叉项 | 实测值 |
| :--- | :---: |
| CID 重叠 | **72 个 = 测试集的 44.7%** |
| 重叠药中标签冲突 | 17（按 TableS1/S3）/ 18（按 Inter-Table4/5），见 §2.1.4 |
| 冲突方向 | **全部是 Y（训练）→ N（测试），零例反向** |
| 重叠药中特征向量不一致 | 3 个（Chloramphenicol 差 10 列、Bupropion 差 38 列、salicylic acid 差 6 列） |

**冲突药名单**：Chloramphenicol, estradiol, Allopurinol, theophylline, Aspirin,
Verapamil, Carbenicillin, domperidone, Ganciclovir, Atropine, ibuprofen, Isoniazid,
progesterone, Propranolol, testosterone, Zidovudine, Doxycycline（17 个）
\+ Bupropion（仅在 Inter 口径下冲突，见 §2.1.4）

### 2.1.1 ⭐ 2×2 列联表 —— 本题的核心证据

72 个重叠药，TableS1（通透性口径）× TableS3（**表名明写 "Based on CNS Activity"**）：

| | S3 = Y | S3 = N |
| :--- | :---: | :---: |
| **S1 = Y** | 39 | **17** |
| **S1 = N** | **0** | 16 |

**`(N, Y)` 格恰好为 0，72/72 = 100% 满足逻辑包含关系。**

这不是噪声的形态。若两表是同一连续量的两次有噪测量，两个非对角格都应有计数。
**一个格恰好为零，是"蕴含"的签名，不是"测量误差"的签名**：
CNS 活性 ⟹ 必然穿过 BBB，故"不穿过却有 CNS 活性"在物理上近乎不可能，实测 0 例。
⇒ S3 的 Y 集，在重叠药上是 S1 的 Y 集的**真子集**。

### 2.1.2 TableS1 内部本身就是三套标准的混合

| 来源文献 | n | Y 占比 |
| :--- | :---: | :---: |
| Li et al., 2005 | 138 | **53%** |
| Abraham et al., 2006 | 53 | **91%** |
| Subramanian & Kitchen, 2003 | 21 | 86% |
| Wang et al., 2015 | 1 | 0% |
| *（对照）*TableS3 / Doniger 2002 | 161 | 47% |

Li 53% 对 Abraham 91%：**同一个名义标签 "BBB+"，正例率差 38 个百分点**。
⇒ 训练集的 213 个标签来自三篇各有自己 cutoff 的论文，表里只记了**引文**，没记 cutoff。

**⚠ 混杂因素（未排除）**：该差距也可能源于各论文**选药范围不同**（Abraham 的集合可能本就
富集 CNS 药），而非 cutoff 不同。故此表只能作为"标签体系异质"的**一致性证据**，
不能单独证明阈值不同。

### 2.1.3 定量数值在仓库中完全不存在（已验证）

| 表 | 全部列 |
| :--- | :--- |
| TableS1 | `CID` \| `Drug Name` \| `BBB +/-` \| `Reference` —— **无任何定量通透性** |
| TableS3 | `ID` \| `Name` \| `BBB +/-` \| `MW` \| `Volume` \| `surface area` \| `%hydrophil` \| `logP` \| `HLB` |

TableS3 的数值列**全是化学描述符，不是通透性测量值**。全仓库 20 个文件按
`logBB / Kp,uu / CSF / brain-plasma ratio / PAMPA / MDCK / Papp` 等词表检索，
**零命中定量通透性字段**（命中项只有模型性能的 AUC、化学描述符 `%hydrophillic surf.`，
以及 `Cerebrospinal fistula` 一类无关副作用词条）。

⇒ **阈值从数据中不可恢复。** 要知道某个 "Y" 指什么，必须回到 Li 2005 / Abraham 2006 /
Subramanian 2003 / Doniger 2002 四篇原文，各自读它的判定标准。这正是用户的论点：
**没有统一标准，而数据集连"用了哪套标准"都只留下一个引文指针。**

### 2.1.4 ✅ 已判明：Bupropion 是仓库真实的内部不一致，不是解析错误

逐单元格回读原始 xlsx（含列字母），四处取值如下：

| 文件 | 单元格 | Bupropion 标签 |
| :--- | :---: | :---: |
| TableS1（附表，训练） | r9 `D` | **Y**（Reference = Li et al., 2005） |
| TableS3（附表，测试） | r8 `C` | **Y** |
| Inter-Table4（特征矩阵，训练） | r10 `CL` | **Y** |
| **Inter-Table5（特征矩阵，测试）** | **r8 `CL`** | **N** ← 与其来源附表 TableS3 相反 |

⇒ **Inter-Table5 与它自己的来源附表 TableS3 矛盾。** 解析无误（列字母明确）。
该药的特征向量在两个文件族间还差 38 列，为全仓库最大，与"该行数据处理出过问题"一致。

**药理学上哪个对**：bupropion 是抗抑郁／戒烟药，**definitively CNS-active**，
按 Doniger 的 CNS 活性判据应为 Y ⇒ **TableS3 的 Y 可辩护，Inter-Table5 的 N 是错的。**

**对出题的影响（重要）**：标准答案一律以**附表 TableS1/TableS3 为标签源**
（有文献出处、可引用），得 17 个冲突与 §2.1.1 的干净 2×2。
**不以 Inter 口径（18）为准** —— 那个多出来的 1 个是仓库录入错误，
不是科学信号。把录入错误写进标准答案会让题目无法辩护。

Inter-Table5 这处矛盾可作为**可选加分子问题**，但不进主答案。

### 2.2 陷阱是什么

#### 2.2.0 根因：BBB 通透性没有统一标准（用户论点，2026-11-05）

用户原话：**"BBB 的通透性本质上就很难说清楚，因为没有统一标准 —— 渗透率 1%、5%、
还是 10%，还是以观察到疗效为准。"**

这是整道题的科学地基。`BBB+/−` 是对一个**连续量**的二分，而二分点有两重任意性：

| 任意性 | 具体表现 |
| :--- | :--- |
| **阈值任意** | 同一个测量（logBB、CSF/plasma 比）取 1% / 5% / 10% 会给出不同的 Y/N |
| **判据任意** | 更根本：到底以**分子是否穿过**为准，还是以**是否观察到临床疗效**为准 |

⇒ `BBB+` 不是一个事实，而是一个**依赖于定义的构造**。任何把它当客观真值的建模都在拟合
一个不自洽的目标。

#### 2.2.1 ⚠ 对我上一轮推断的修正

上一轮我写"零例反向说明两套标签定义不同"。**用户的论点暴露出这个推理是不充分的**：

| 竞争解释 | 预测的冲突方向 |
| :--- | :---: |
| ① 判据不同（穿过 vs 临床疗效） | 单向 Y→N |
| ② 判据相同，Doniger 的阈值更严 | **同样是单向 Y→N** |

**两种解释预测同一个签名，故单看方向无法区分。** 我原来的推断是欠定的。

**现在支持"判据不同"的证据另有两条，都不依赖方向**：

1. **文献证据（直接）**：TableS3 的表名就写着 **"Independent Test Data Set of Drugs Based on
   CNS Activity"**。判据写在表名里。
2. **结构证据（强）**：`(N, Y)` 格为 **0/72**。若只是阈值更严，被翻转的药应与"是否 CNS 活性"
   无关，不必然产生恰好为零的格；而"CNS 活性 ⟹ 穿过"是一条**逻辑蕴含**，它预测这一格
   **必须**为零。实测 100% 吻合。

**仍未做的核实**：四篇原文各自的判定标准未逐篇回查，故"S1 系 = 穿过口径"这半句仍是推断。
定稿前须读 Li 2005 / Abraham 2006 / Subramanian 2003 / Doniger 2002 的 Methods。

#### 2.2.2 陷阱的机制 —— 同一份数据，两个指标，只有一个能看出答案

| | 做法 | 得到 | 能否看出真相 |
| :--- | :--- | :---: | :---: |
| **非专家** | 数"有多少个不一致" | 17（≈24% 分歧率） | **不能** —— 看起来就是脏数据 |
| **专家** | 建完整 2×2 列联表 | 39 / 17 / **0** / 16 | **能** —— 那个零说明一切 |

**这是本题最漂亮的地方：两个指标算自同一份数据，成本几乎相同，但只有后者能揭示答案。**
分歧计数把"结构"压成了"比例"，而答案恰恰只存在于结构里。

**非专家会怎么做**（四条都是错的）：合并两表取并集；多数投票；信"更新的那一篇"；
把 17 个冲突当噪声丢掉。
**专家会怎么做**：先建 2×2，看到结构零，判定两表是**蕴含关系**而非**重复测量**，
于是保留来源作为分层变量，并认定那 17 个药是全数据集**信息量最大**的样本 ——
它们恰好就是"穿过了 BBB 但不是 CNS 治疗药"的那一类。

### 2.3 Decision table（§2A 要求，定稿前须经用户确认）

| 项 | 本题的定义 |
| :--- | :--- |
| **研究问题** | 这两张表的 BBB 标签可否视为同一真值并合并用于训练/评估？ |
| **观察/证据单位** | 一个 (CID, 来源表) 对 —— **不是** 一个药物。同一药物在两表中是两条独立证据。 |
| **保留与去重规则** | 按 CID 对齐，但**不得合并、不得择一保留**；冲突条目必须双方都保留并可追溯到来源表、来源文献、原始标签值。 |
| **判定指标** | ① 重叠 CID 数；② **完整 2×2 四格计数（核心，答案就在 `(N,Y)=0` 这一格）**；③ 冲突数（仅作对照，单独看它不足以答对）；④ TableS1 内各来源的 Y 占比；⑤ 特征向量不一致数。 |
| **允许的结论** | 可结论"两表是蕴含关系而非重复测量，不可合并、不可多数投票、不可丢弃冲突"；可结论"阈值从数据不可恢复"。**不可**结论"某一方标签错误"；**不可**在未回查四篇原文 Methods 前断言各自的确切判定标准；**不可**仅凭 §2.1.2 的来源间差异断言阈值不同（选药范围混杂未排除）。 |
| **过滤影响审计** | 须输出：按常规做法丢弃冲突药后，训练集与测试集各损失多少条、Y/N 比如何变化、以及被丢掉的恰是哪 17 个药（须给药名，因为它们正是信息量最大的样本）。 |

### 2.3.1 Grader 设计（对上 LatchBio 公开格式）

| 字段 | 类型 | 标准答案 | Grader |
| :--- | :---: | :---: | :---: |
| `n_overlap` | int | 72 | NumericTolerance, tol=0 |
| `cell_YY` / `cell_YN` / `cell_NY` / `cell_NN` | int | 39 / 17 / 0 / 16 | NumericTolerance, tol=0 |
| `relation` | choice | `implication`（而非 `noisy_replicate`） | MultipleChoice |
| `correct_handling` | choice | `stratify_by_source`（而非 merge / vote / drop） | MultipleChoice |

**为什么这个 grader 稳健**：四格计数是**确定性整数计数**，不涉及 normalization、随机种子或
阈值选择，因此 tol=0 也完全可复现。这正好避开了 LatchBio 自己在
`Human Verification of SpatialBench` 里报告的最大问题 —— 他们有一道题的相关系数因 normalization
不同从 0.62 变成 0.85，容差定得太窄导致专家答案被误判为错。
**本题从源头上不存在那个自由度。**

### 2.3.2 ⭐ 四篇来源论文的判定标准（2026-11-05 文献核实，§4 第 2 步已执行）

| 来源 | 判据 | 阈值 | 核实程度 |
| :--- | :--- | :---: | :---: |
| **Doniger et al., 2002**（= TableS3） | **CNS 活性**（active / inactive in CNS），**不是** logBB | **无 logBB 阈值** | **已证实** |
| Subramanian & Kitchen, 2003 | logBB（log 脑/血浆浓度比） | logBB ≥ 0.3 / < −1 | **仅二手来源，原文未核** |
| Abraham et al., 2006 | logBB；**假设实验误差 ≈ 0.3 log 单位**，logBB 范围 −2.2 ~ 1.6 | ≥ 0.3 易穿过；−1 ~ 0.3 仍可穿过；< −1 难入脑 | 部分；**具体是哪一篇未定**（见下） |
| **Li et al., 2005** | 分类型 BBB+/BBB− | **未查到** | **未解决，PubMed 被 cookie 墙拦住** |

**TableS3 = Doniger 数据集，已由三条独立证据证实**：
1. 表名：`Independent Test Data Set of Drugs Based on CNS Activity`
2. **列名与 Doniger 的 9 个描述符逐一对应**：TableS3 给出 MW、Volume、surface area、
   %hydrophilic surf.、logP、HLB —— 正是 Doniger 九项中的 6 项，顺序一致
3. Doniger 原文自述判据为 CNS 活性（179 BBB+ / 145 BBB−，共 324 化合物）

⇒ **§2.2.1 里"S3 = CNS 活性口径"这半句，已从推断升级为事实。**

### 2.3.3 ⭐⭐ 用户论点的定量证据 —— 文献里同时在用五个不同阈值

用户原话问"渗透率 1%、5%、还是 10%"。文献实际在用的 logBB 阈值换算成脑/血浆浓度比：

| 文献采用的 logBB 阈值 | 脑/血浆浓度比 | 出处 |
| :--- | :---: | :---: |
| ≥ +0.3 = BBB+ | **200%** | Abraham 系、Kunwittaya 等 |
| 0 ~ −1（"通常最优区间"） | 100% ~ 10% | 多篇综述 |
| ≤ −0.3 = BBB− | **50%** | Martins / Wichmann 系 |
| −0.52（"更合适的 BBB− 界"） | **30%** | 后续修订意见 |
| < −1 = BBB− | **10%** | Kunwittaya 等 |

**⇒ 同一个名义标签 `BBB+`，界线在 10% 到 200% 之间横跨 20 倍。** 用户的直觉方向正确，
实际离散程度比"1%/5%/10%"更大。

**三条更锋利的推论（均有文献支持）**：

1. **阈值等于实验误差。** BBB+ 的界线取 logBB ≥ 0.3，而这个 0.3 正是 Abraham 假设的
   **实验误差本身**（≈0.3 log 单位）。⇒ 单次 logBB 测量**无法可靠地把化合物定位在它自己的
   分类边界的哪一侧**。部分数据集（Martins/Wichmann）因此直接**剔除 −0.3 ~ +0.3 之间的化合物**，
   理由是"缺乏明确的 permeable/non-permeable 性格"。
2. **模型误差已超过这个边距。** B3DB 中 1058 个化合物有 122 个 `|logBB| < 0.1`，
   而最优测试 RMSE ≈ **0.37 log 单位 > 0.3 的边距**。⇒ 二分边界比实验误差和模型误差都更细。
3. **连底层连续量本身都不固定。** logBB 测量**高度依赖采样时间** ——
   所以它不是药物的一个固定属性。

**最关键的一条（直接支撑 Task 1 的立意）**：文献已明确记载，
**"logBB 阈值变化会改变预测准确率，故基于 logBB 的二分 BBB 模型易引入阈值设定偏倚，
这正是采用 CNS-活性数据集的动因 —— 后者不需要设阈值，从而避开该偏倚。"**

⇒ **TableS1（logBB 族）与 TableS3（CNS 活性族）不是同一事物的两次测量，
而是这个领域一个已知方法学困境的两只角。** 这已是文献记载，不再是我的推断。

### 2.4 ⚠ 必须由用户拍板的一点

Task 1 的标准答案会指出**该公开数据集的"independent test set"与训练集有 44.7% 重叠**。
这是用户本人论文的附带数据。

- **我未读该论文 Methods**，因此**不能断言论文的建模确实存在 train/test leakage** ——
  论文可能在训练时已剔除重叠药，附表只是原始数据罗列。此点**未验证**。
- **建议的规避写法**：题目不问"这篇论文对不对"，只给两张表、问"能否合并、应如何正确评估"。
  这样标准答案不依赖论文有无缺陷，用户也无需评价自己的工作。
- 若用户愿意直接用 leakage 框架出题，需先一起读论文 Methods 确认事实。

---

## §3 Task 2 — glioma 药物治疗（仅设计方向，未开工）

### 3.1 用户的科学主张（出题的立意来源，用户原话转述）

glioma 化疗首先须突破 BBB；但晚期 BBB 已自我破坏；低氧 stress 造成癌细胞生长停滞或死亡；
而**生长停滞本身又增强放疗与化疗抗性**。

### 3.2 ✅ 两点收紧已回查文献核实（2026-11-05，§4 第 4 步）

| 需收紧处 | 文献结论 | 核实 |
| :--- | :--- | :---: |
| "晚期 BBB 才破坏" | **强确认且更有利**：GBM 的屏障是**空间异质的马赛克** —— 核心血管通透性升高、屏障破坏，而**浸润边缘保持基本完整的 BBB 特性，从而保护浸润瘤细胞免受药物作用**；且明确记载**复发正源于边缘** | ✅ |
| "生长停滞→抗性"非唯一机制 | 确认：低氧致放疗抗性经**氧效应**（氧把辐射所致 DNA 自由基固定为永久损伤）；glioma 实测 OER ≈ **1.23–1.54**，无氧极限下约需 **3 倍剂量**；TMZ 抗性另涉 MGMT 与 MMR | ✅ |

**额外收获（文献给的比原设计更好的陷阱）**：
- **增强程度不是药物暴露的有效代理**：屏障状态是马赛克，MRI 上的"渗漏"是"有效药物暴露"的
  **拙劣代理**。
- **bevacizumab 的悖论**：抗 VEGF **通过恢复屏障完整性来减少增强** —— 影像学上像是好转，
  药理学上却是新障碍；此后复发**特征性地不增强、且远离原核心**。
- **低氧维持静止态的干细胞样群体**，该静止组分**直接关联治疗抗性与复发**。

### 3.3 ⚠️ Task 2 的 v1/v2 设计（**已被 §3.6 的 v3 取代，保留作对照**）

> **本节描述的不是当前有效版。** 当前有效版是 §3.6 的
> [evals/glioma_margin_rate_limit_v4.json](../evals/glioma_margin_rate_limit_v4.json)（§3.7）。
> 本节所述设计因"不需要 glioma 知识即可全对"于 2026-11-05 废弃，归档在
> [superseded_for_comparison/evals/v2_20261105_elimination_solvable/](../superseded_for_comparison/evals/v2_20261105_elimination_solvable/)。
> 保留本节是因为 v3 的每一处改动都是针对本节的某个具体缺陷，删掉就看不出改了什么。

**文件**：`tmp/glioma_margin_rate_limit_v1.json`（**已废弃**，2026-11-05 按用户指示从 `evals/` 移入 `tmp/`）

**题干**：复发性 GBM，给三个药及其**本数据集中的 BBB 标签**，加四条临床事实
（核心屏障破坏、边缘屏障完整、复发源于边缘、边缘有低氧静止组分）。
先问现状，**再假设一项递送技术完全消除边缘屏障**，重新评估。

| 药 | 数据集标签 | 机制 |
| :--- | :---: | :--- |
| lomustine | **Y** | 亚硝基脲，烷化并交联 DNA |
| methotrexate | N | 抗叶酸，**S 期特异** |
| doxorubicin | N | 拓扑异构酶 II 毒物，**需复制**才能转为致死双链断裂 |

**六个答案字段全为类别** ⇒ MultipleChoiceGrader，无容差需校准。

**核心陷阱**：问"递送解决后还剩什么限制"，要求选**三药共有性质**所致的那一项。
答案 `target_cell_quiescence` —— **三药的细胞毒性都依赖 DNA 复制**，
静止细胞不把任何一种药-靶相互作用转导为死亡，药到多少都一样。

⇒ **用户那条因果链的悖论在此落地**：低氧 stress 造成的生长停滞，
表面上像肿瘤在承受代谢压力，**实际正是使其难治的那个状态**。

**最强干扰项 = `mgmt_mediated_repair`**：MGMT 甲基化是 GBM 被引用最多的抗性机制，
模型按"GBM 耐药"做模式匹配会首先取它。它在此错在两点：
① 只针对 O6-烷基鸟嘌呤，**解释不了 methotrexate 与 doxorubicin 的失败**；
② 题目明问"三药共有性质"。`mismatch_repair_deficiency` 是同类陷阱。

**为什么不用 temozolomide**：它是标准治疗、也是复制依赖性最清楚的例子，
但**本数据集里没有它**（已检索确认）⇒ 用它就得断言一个数据不承载的 BBB 标签。
这与 Task 1 "标签必须可追溯到来源"的立意自相矛盾，故舍弃。

### 3.3 候选题设计（倾向方案）

**题干**：给定一组候选药的 BBB 通透性标签与体外抗增殖活性，问哪一个最不值得推进到
复发性 glioma，并给出理由类别。

**陷阱（三层，答对需全部穿过）**：
1. 最直觉的答案是"BBB 不通透的那个"——但肿瘤核心 BBB 已破坏，强化区可达药。
2. 第二层是"体外最强的那个最好"——但体外活性读数受细胞系倍增时间影响。
3. **真正的答案**：机制依赖细胞增殖的药，对低氧区生长停滞的细胞无效，
   而那群细胞正是复发来源 —— **BBB 通透且体外最强的药，可能恰恰最不值得推进。**

**答案形式**：类别（选项）+ 理由类别，对应 MultipleChoiceGrader。
这比数值容差稳健，不同合理流程不会给出不同答案。

### 3.4 ⭐ 用户的"无统一标准"论点同时解掉了 Task 2 的数据源阻塞

原本 Task 2 卡在"体外活性数据没有用户熟悉的来源"。§2.2.0 的论点提供了一条更好的出路：

**把 Task 2 建成 Task 1 的生物学续章，共用同一份数据集，不引入第二个数据源。**

| | Task 1 | Task 2 |
| :--- | :---: | :---: |
| 层面 | 数据层 | 生物学层 |
| 问的是 | `BBB+` 这个标签**依赖于定义** | 在 glioma 里，这个标签**根本不是限速步骤** |
| 答对的关键 | 看出 2×2 的结构零 | 看出问题被问错了 |

**两道题合起来是一个连贯的主张**：BBB 通透性既不是客观事实（Task 1），
在 glioma 这个具体场景里也不是决定成败的那个变量（Task 2）。

**Task 2 的三层陷阱可以完全不依赖药敏数值**：
1. 直觉答案"选 BBB 通透的药" —— 但在 glioma，**肿瘤核心的 BBB 已破坏**，强化区可达药。
2. 进一层"那就无所谓了" —— 但**浸润边缘的 BBB 相对完整**，而复发正源于边缘的残留细胞。
3. **真正的答案**：边缘残留细胞处于低氧、生长停滞状态，**而生长停滞本身即赋予放化疗抗性**。
   ⇒ 限速步骤不是"药能不能进去"，而是"进去之后那群细胞是否处于药物机制所需的状态"。
   机制依赖增殖的药，即使 BBB 通透且体外最强，也打不到复发来源。

**答案形式**：类别选择（限速步骤是哪一个）+ 理由类别 ⇒ MultipleChoiceGrader，无容差问题。

### 3.5 Task 2 仍未解决项

| 阻塞 | 状态 |
| :--- | :--- |
| §3.2 两处收紧未回查文献 | **✅ 已解决**（§3.2，两条均确认，且另得三条更好的陷阱） |
| 候选药清单未定 | **✅ 已解决**（lomustine / methotrexate / doxorubicin，标签均经本数据集实测） |
| 题面如何不泄题 | **✅ 已修（v1→v2，用户 2026-11-05 选定方案 A）**。v1 的事实 #4 写了 `growth-arrested`、字段名含 `arrest`，与正确选项 `target_cell_quiescence` 近义 ⇒ **可仅靠词形匹配选中，不需任何生物学**。v2 把事实 #4 减为"存在低氧 niche"、字段改名 `hypoxia_effect_on_cytotoxic_sensitivity`。改后 `arrest` 在全题面消失，`quiescen` 只剩选项列表 1 处（必须保留）。作答所需前提仍全在题面（低氧 niche + 各药机制），删掉的只是本该由作答者完成的推理。v1 已归档并可并排 diff：[superseded_for_comparison/evals/v1_20261105_lexical_leak_growth_arrested/](../superseded_for_comparison/evals/v1_20261105_lexical_leak_growth_arrested/) |
| **题目整体太简单（2026-11-05 自审 + 用户追问）** | **✅ 已重做为 v3，见 §3.6**。v2 六字段中四个白送/通识，决定性字段可纯排除法解出 |
| 强模型可能直接答对 | ⬜ **未实测，这是最大的未知**。v3 的方向反转字段也未实测 |
| 无独立专家预先作答 | ⬜ 未做。LatchBio 自己的人工核验表明，**未经独立专家从题面出发作答的题不可信**（159 题首轮仅 94 题通过，部分失败源于规格不清而非能力不足） |

### 3.6 ⚠️ Task 2 的 v3 设计（**已被 §3.7 的 v4 取代，保留作对照**）

> **本节描述的不是当前有效版。** 当前有效版是 §3.7 的
> [evals/glioma_margin_rate_limit_v4.json](../evals/glioma_margin_rate_limit_v4.json)。
> v3 因**题面在前提里陈述了该考的疾病生物学**于 2026-11-05 废弃，归档在
> [v3_20261105_premises_stated_the_biology/](../superseded_for_comparison/evals/v3_20261105_premises_stated_the_biology/)。
> 本节保留，因为 v4 的每处改动都针对本节某个具体缺陷。

**文件**：`tmp/glioma_margin_rate_limit_v3_OBSOLETE.json`（已废弃，按用户指示移入 `tmp/`）

#### 3.6.1 为什么推翻 v2 —— 两条独立理由

**理由 A（我自审出的结构缺陷）**：v2 六个字段里四个是题面白送或教科书通识；
决定性字段 `remaining_primary_limit` 的四个干扰项**全是药物特异**，正确选项是唯一药物无关的，
⇒ **纯排除法可解，不需要任何 glioma 知识**。而题面里我自己写的那句
`Choose the factor that applies to all three agents by virtue of a property they share`
**主动给这条排除路径立了路标**。

**理由 B（用户 2026-11-05 指出，更根本）**：
> 「BBB 在脑瘤中不是特别重要，因为 glioma 绝大多数伴随 BBB 破坏，可以作为诱骗项。」

v2 把"穿不过 BBB"当成了真实障碍。但增强区 BBB 早已自我破坏 ⇒ 药物的完整-BBB 通透性
在那里不决定暴露；而在屏障完整的边缘区，暴露又不决定杀伤（细胞不增殖）。
⇒ **`BBB +/-` 是数据集里最显眼、却不起作用的那个变量，它应该是诱骗项。**
v2 反而把它当限速步骤，等于主动放弃了最好的陷阱。

**由此得到的一般判据（可复用）**：
**凡是正确答案与文献主流叙事同向的陷阱都不可能难** —— 模型检索到主流叙事就对了。
难的陷阱必须满足：**主流叙事给出的答案是错的。** v2 的低氧→静止→复制依赖链在 glioma
综述里是标准叙事，所以它注定简单。

#### 3.6.2 v3 的五个字段（6 → 5，每个都独立鉴别）

| 字段 | 标准答案 | 它拦住什么 |
| :--- | :---: | :--- |
| `primary_determinant_of_margin_response` | `proliferative_state_of_the_target_cells` | **BBB 诱骗项**；五选项全为药物无关，排除法失效 |
| `agent_selection_by_bbb_label` | `lomustine_preferred_but_this_does_not_predict_response` | 把表读对仍可能答错 —— 问的是读法有无预测力 |
| `antivegf_effect_on_margin_exposure` | `decreases_via_barrier_restoration` | **方向反转**：vascular normalization 框架预测相反 |
| `antivegf_effect_on_margin_cytotoxic_sensitivity` | `decreases` | **方向反转**：联用的已发表理由预测相反 |
| `relapse_pattern_after_antivegf` | `non_enhancing_and_distant` | **自洽性检验**：方向答反者在此会自相矛盾 |

#### 3.6.3 逐条改动对应哪个缺陷

| 改动 | 针对 |
| :--- | :--- |
| 删掉 `...by virtue of a property they share` 提示句 | 理由 A |
| 决定性字段五选项全改为药物无关（屏障通透性 / 增殖状态 / 血管外间质转运 / DNA 修复能力 / 免疫微环境） | 理由 A |
| 删两个白送字段、合并两个重复 hypoxia 字段 | 理由 A |
| 新增两个方向反转字段（持续抗 VEGF 对暴露、对内在敏感性） | 主流叙事判据 |
| 新增 `agent_selection_by_bbb_label` | 理由 B |
| 新增 `relapse_pattern_after_antivegf` 作自洽性检验 | 主流叙事判据 |

#### 3.6.4 ⚠ 为了可辩护性而必须加的限定

抗 VEGF 存在**短暂的血管正常化窗口**，窗口内氧合与递送确实改善。
⇒ 题面前言与每个抗 VEGF 字段描述都**明确限定为"超出该窗口之后的持续治疗"**。
**不加这个限定，`increases_via_reduced_interstitial_pressure` 就是站得住的，题目会变成"有歧义"而不是"难"。**
这个限定虽然提示了窗口的存在，但不提示方向，是必要的交换。

#### 3.6.5 v3 的已知缺口（已写入 eval 的 `known_gaps`，不隐藏）

| 缺口 | 影响 |
| :--- | :--- |
| 持续抗 VEGF 下低氧加深，更多依赖临床前模型 | 答案方向站得住，但**是推断不是测量**；临床侧证据是复发表型，不是边缘区氧电极 |
| 无来源直接测量这三药对不增殖 glioma 细胞的 sparing | 故措辞落在"复制依赖性"，不断言 sparing 比例 |
| Xenopus 提取物提出烷化剂有不依赖复制的成分 | **减弱但不反转**结论 |

#### 3.6.6 v3 仍未解决

- **未在任何模型上实测。** v2 的缺陷是我自审出来的、不是实测出来的；v3 是否真更难同样需实测。
- **方向反转字段的残余风险**：若某模型恰好读过 GBM 里抗 VEGF 屏障恢复那一支文献，会直接答对。
  此时难度来自知识覆盖而非推理，仍需进一步加深。
- ~~`evals/glioma_margin_rate_limit_v1.json`（已废弃）仍留在 `evals/` 中~~
  **✅ 已解决（用户 2026-11-05 指示"过时的文件转入 tmp"）**：已 `mv` 至 `tmp/`，移动前经 `cmp`
  核对归档件逐字节相同。⇒ **`evals/` 内现只有两个交付 JSON + README，无非交付文件。**
  见 [evals/README.md](../evals/README.md) 与 [tmp/README.md](../tmp/README.md)。

### 3.7 ⭐⭐⭐ Task 2 当前有效设计 = v4（2026-11-05，删掉前提里的全部疾病生物学）

**文件**：[evals/glioma_margin_rate_limit_v4.json](../evals/glioma_margin_rate_limit_v4.json)
（英文渲染 [tmp/glioma_margin_rate_limit_v4_1105.md](../tmp/glioma_margin_rate_limit_v4_1105.md)、
中文审阅译文 [tmp/glioma_margin_rate_limit_v4_zh_1105.md](../tmp/glioma_margin_rate_limit_v4_zh_1105.md)）

#### 3.7.1 为什么推翻 v3（用户 2026-11-05 提问）

> 「为什么要告诉 AI 这些事实？ 该肿瘤具有以下事实：对比增强的肿瘤主体 BBB 已被破坏，血管通透性升高。」

v3 题面有一个四条事实的清单，**四条全是 GBM 一般生物学、不是这个患者的情况**：

| v3 的"事实" | 性质 | 判定 |
| :--- | :---: | :---: |
| 1. 增强主体 BBB 已破坏、通透性升高 | GBM 通识 | **删（致命的那条）** |
| 2. 非增强浸润边缘屏障基本完整 | GBM 通识 | **删** |
| 3. 复发源自边缘浸润细胞 | GBM 通识 | **删** |
| 4. 边缘区含低氧 niche | GBM 通识 | **删** |

**第 1 条自己拆了自己的陷阱**：v3 的立意是"BBB 标签是诱骗项"，而诱骗项有效**全靠模型自己不知道标签不管用**；
写上第 1 条 = 直接告诉它标签不管用。

v3 还逐药写了机制（`S-phase specific` / `during replication` / `alkylates and crosslinks DNA`），
那是标准药理学，等于供出了**复制依赖性**这一半推理。

#### 3.7.2 ⚠⚠ 这是同一个错误的第三次 —— 根因与可复用判据

| 版本 | 泄题形式 | 层级 |
| :--- | :--- | :---: |
| v1 | 前提写 `growth-arrested`，字段名含 `arrest` | **词形** |
| v2 | "选择因三药共有属性而成立的因素" ⇒ 给排除法立路标 | **解法** |
| v3 | 四条事实陈述疾病生物学 ⇒ 给出推理所需全部前提 | **知识** |

**根因**：我把 LatchBio `adding_evals.md` 的"充分指定"误解成"前提喂到嘴里"。正确界线：

> **患者层面、因病例而异、无法先验推知的 ⇒ 必须写明。**
> **疾病层面、对所有 GBM 都成立的 ⇒ 一律不写，那正是考题。**

⇒ 该原则已写成 v4 的 `notes.specification_principle`，是**题目的显式规格**而非我的记性。
配合 §3.6.1 那条判据（**正确答案与文献主流叙事同向的陷阱不可能难**），两条合起来是本项目的出题准则。

#### 3.7.3 v4 的改动

| 改动 | 效果 |
| :--- | :--- |
| 整块删除四条事实清单 | 屏障状态、复发起源、低氧全部交还模型 |
| 删除三药机制描述，只留药名 + 数据集标签 | 复制依赖性交还模型 |
| 两个抗 VEGF 字段选项缩为**裸方向** | 原选项 `unchanged_because_margin_barrier_already_intact` **本身就陈述了边缘屏障完整**，是第五处泄题 |
| 新增 `notes.specification_principle` | 防第四次 |
| `known_gaps` 新增一条 | 诚实记录 v4 引入的新失败模式（见 §3.7.5） |

**题面长度：2799 → 2186 字符。** 五个字段与标准答案不变（两个抗 VEGF 字段的答案由
`decreases_via_barrier_restoration` 简化为 `decreases`）。

#### 3.7.4 实测复查：题面已无疾病生物学残留

```bash
cd /home/gao/projects_2026Q4/48_Saul_LatchBio
python3 -I -c "
import json,re
t=json.load(open('evals/glioma_margin_rate_limit_v4.json'))['task']
for w in ['disrupted','intact','hypox','niche','recurrence arises','S-phase','replication',
          'alkylat','crosslink','topoisomerase','antifolate','nitrosourea','mosaic']:
    print('%-18s %d'%(w,len(re.findall(w,t,re.I))))
"
```

2026-11-05 实测：以上**全部为 0**。`enhanc` 2 处命中均在 `relapse_pattern_after_antivegf`
的选项名内（`enhancing_and_local` / `non_enhancing_and_distant`）—— **是答案选项而非前提，必须保留。**

#### 3.7.5 ⚠ v4 新引入的代价（已写入 `known_gaps`，不隐藏）

不再陈述药物机制 ⇒ **不知道 methotrexate 是 S 期特异、或不知道 doxorubicin 致死需要复制的模型，
会因为缺标准药理学而失败，而不是因为缺预期的洞察。** 选这三个药的部分理由正是这套药理学并不冷僻，
但这个失败模式真实存在，读实测结果时必须把它与"真的没想到"区分开。

另：为可辩护性保留的"超出正常化窗口之后"那句限定仍在（§3.6.4），它暴露窗口的存在、但不暴露方向。

#### 3.7.6 v4 仍未解决

- **未在任何模型上实测。** v1→v2→v3→v4 四代缺陷**全部**来自自审或用户追问，**没有一条来自实测** ——
  这意味着当前对 v4 难度的全部主张都仍是判断而非证据。**这是现在最该做的一件事。**
- 方向反转字段的残余风险：若模型恰好读过 GBM 抗 VEGF 屏障恢复那一支文献，会直接答对；
  此时难度来自知识覆盖而非推理。
- 废弃版本均已移入 `tmp/`，`evals/` 内只剩两个交付 JSON。目录内哪个有效见
  [evals/README.md](../evals/README.md) 与 [tmp/README.md](../tmp/README.md)。

---

## §4 下一步（待用户决策后执行）

| # | 动作 | 前置条件 |
| :---: | :--- | :---: |
| # | 动作 | 状态 |
| :---: | :--- | :---: |
| 1 | 判明 Bupropion 的仓库内部不一致 | **✅ 已完成**（§2.1.4，确为仓库错误，标准答案改以附表为准） |
| 2 | 读四篇来源论文的判定标准 | **🟡 3/4 完成**（§2.3.2；Doniger 已证实，Li 2005 未查到） |
| 3 | 确认仓库 license | **✅ 已完成**（§1.1，无 license，用户加一个即可解除） |
| 4 | 回查 §3.2 两处 glioma 文献收紧 | **✅ 已完成**（§3.2，两条均确认 + 三条额外陷阱） |
| 5 | 决定 Task 1 框架 | **✅ 用户已拍板：蕴含关系框架**（2026-11-05） |
| 6 | 写两道题的题面 + grader + 标准答案复现脚本 | **✅ 已完成**（见 §4.2） |
| 7 | 在 2–3 个模型上各跑数次实测，留原始输出 | ⬜ **下一步，最大未知** |
| 8 | 向 Saul 确认 90 min 与可否预先准备 | ⬜ **用户发送**（只准备不发送） |

### 4.2 已产出的交付物

| 文件 | 作用 |
| :--- | :--- |
| [evals/bbb_label_implication_v1.json](../evals/bbb_label_implication_v1.json) | Task 1，NumericTolerance，7 字段，全部 tol=0 |
| [evals/glioma_margin_rate_limit_v4.json](../evals/glioma_margin_rate_limit_v4.json) | **Task 2 当前有效版**，MultipleChoice，5 字段（v1–v3 已废弃，见 §3.7） |
| [evals/README.md](../evals/README.md)、[tmp/README.md](../tmp/README.md) | 声明哪个是当前有效版（目录内有多版本时的必备配套） |
| **[Submission/submission_both.txt](../Submission/submission_both.txt)** | ⭐ **限时作答时直接用这份**：两题合一 390 行，答案紧跟题目 |
| [Submission/submission_task1.txt](../Submission/submission_task1.txt)、[Submission/submission_task2.txt](../Submission/submission_task2.txt) | 单题版，若输入框按题分开则用这两份 |
| [scripts/1_task1_ground_truth.py](../scripts/1_task1_ground_truth.py) | Task 1 标准答案复现，**仅标准库**，确定性 |
| [scripts/task1_ground_truth.json](../scripts/task1_ground_truth.json) | 脚本产出的机器可读答案 |
| docs/（本文）、Email.txt、LatchBio.txt | 计划与原始邮件 |

**两个 eval 已过 JSON 语法与字段校验**：`id` / `task` / `data_node` / `grader` / `notes`
五个必需键齐全，`notes` 下 `platform` / `grader` / `what_this_tests` / `key_metrics` /
`challenges` 五个推荐子键齐全。

### 4.3 ⚠ 交付前仍必须做的两件事

| # | 事项 | 为什么不能省 |
| :---: | :--- | :--- |
| 1 | **双向验证**：正确答案必须通过、错误答案必须不通过 | LatchBio `adding_evals.md` 明文要求；"从未用错误答案测过的 eval 不可靠" |
| 2 | **模型实测** | 若强模型直接答对，陷阱即失效，须换陷阱。**这是目前最大的未知** |

**Task 1 的 grader 组合需注意**：LatchBio 格式为"一个 eval 一个 grader type"，
而 Task 1 含 5 个整数字段 + 2 个类别字段。已在 `notes.grader` 里写明：
若其 harness 要求单一 type，**仅 5 个整数字段即足以区分对错**，
两个类别字段可拆成配套的 MultipleChoice eval。

### 4.1 三项仍未闭合的文献缺口

| 缺口 | 影响 | 谁能最快解决 |
| :--- | :--- | :---: |
| **Li et al., 2005 的判据未查到** | 它占训练集 138/213 = **65%**，是最大单一来源。其判据未明 ⇒ 不能断言"TableS1 整体 = logBB 口径" | PubMed 被 cookie 墙拦住；**用户当年选了这篇，可能直接知道** |
| Subramanian & Kitchen 的阈值仅二手 | §2.3.3 的 ≥0.3 / <−1 归于该文，但未见原文 | 需原文 PDF |
| "Abraham et al., 2006" 是哪一篇未定 | 可能是 Abraham 2004 *Eur J Med Chem*，或 Platts/Abraham 等的 LFER 大数据集研究 | **用户当年引的，最清楚** |

**这三项不阻塞出题**：Task 1 的标准答案建立在 2×2 的结构零上，而该结构零已由
TableS3 = CNS 活性口径（已证实）充分解释。上述缺口影响的只是"TableS1 那一侧的口径
能否逐篇说清"，属于题目注释的完备性，不影响答案正确性。

**§8 自审纪律**：每道题须验证"正确答案通过 / 错误答案不通过"（LatchBio 公开文档
`adding_evals.md` 的明文要求），且标准答案须记录软件版本、随机种子、参数。
