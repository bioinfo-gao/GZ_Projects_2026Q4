# glioma_margin_rate_limit_v1（中文审阅版）

> ⚠ **本文是译文，不是真身。** 内容真身在 `evals/glioma_margin_rate_limit_v1.json`（英文）。
> 要改内容请改那个 JSON，然后重译本文。交付给 LatchBio 的是英文 JSON，本文不交付。
> 翻译日期：2026-11-05 ｜ 对应版本：v2（已修字面泄题，v1 归档在
> `superseded_for_comparison/evals/v1_20261105_lexical_leak_growth_arrested/`）

**数据：** https://github.com/bioinformatics-gao/CASE-BBB-prediction-Data

## Grader — MultipleChoice

| 字段 | 标准答案 | 容差 |
| :--- | :---: | :---: |
| `agent_reaching_margin_unaided` | lomustine | 精确 |
| `enhancement_as_exposure_proxy` | invalid | 精确 |
| `delivery_alone_sufficient` | no | 精确 |
| `remaining_primary_limit` | target_cell_quiescence | 精确 |
| `hypoxia_effect_on_cytotoxic_sensitivity` | decreases | 精确 |
| `hypoxia_effect_on_radiotherapy` | reduces_efficacy | 精确 |

## 题面（交给模型的原文，此处为中译）

一位患者患复发性胶质母细胞瘤（glioblastoma）。你要评估三种细胞毒药物，它们的血脑屏障（blood-brain barrier, BBB）通透性标签记录在随附数据集中（`1) Supplementary-data/TableS1-Drug-for-Train-and-Validation.xlsx`，列 `BBB +/-`）：

* **lomustine** —— 数据集标签 `Y`；亚硝基脲类，烷基化并交联 DNA
* **methotrexate** —— 数据集标签 `N`；抗叶酸药，阻断胸苷酸合成，S 期特异
* **doxorubicin** —— 数据集标签 `N`；拓扑异构酶 II 毒物，在复制过程中把酶-DNA 复合物转化为双链断裂

该肿瘤具有以下临床与组织学事实：

1. 对比增强的肿瘤核心区 BBB 已被破坏，血管通透性升高。
2. 非增强的浸润性边缘区 BBB 性质基本完整。
3. 复发源自边缘区的浸润细胞，而非核心区。
4. 边缘区含有一个低氧（hypoxic）niche。

请回答下列问题。然后假设有一种假想的递送技术，可在浸润性边缘区完全消除 BBB，使三种药物到达边缘组织的浓度与到达核心区相同。重新评估。

把答案以 JSON 写入 `eval_answer.json`，字段恰好为：

* `agent_reaching_margin_unaided`（字符串）—— `lomustine`、`methotrexate`、`doxorubicin`、`none` 之一。基于数据集标签，三者中哪一个预期能在没有该递送技术的情况下到达浸润性边缘区。
* `enhancement_as_exposure_proxy`（字符串）—— `valid` 或 `invalid`。MRI 上的对比增强范围，是否是驱动复发的那群细胞所受**有效药物暴露**的有效替代指标。
* `delivery_alone_sufficient`（字符串）—— `yes` 或 `no`。在该递送技术就位、渗透不再是限制因素后，这些药物是否预期能根除驱动复发的细胞群。
* `remaining_primary_limit`（字符串）—— `target_cell_quiescence`、`mgmt_mediated_repair`、`drug_efflux_bcrp`、`insufficient_dose`、`mismatch_repair_deficiency` 之一。渗透问题解决后，对三种药物**共同**仍起主要限制作用的那个因素。请选择那个因为三者共有的某一属性而对三者都成立的因素。
* `hypoxia_effect_on_cytotoxic_sensitivity`（字符串）—— `increases`、`decreases`、`no_effect` 之一。边缘区的低氧状态对那些细胞对这三种药物敏感性的影响。
* `hypoxia_effect_on_radiotherapy`（字符串）—— `increases_efficacy`、`reduces_efficacy`、`no_effect` 之一。同一低氧状态对同期放射治疗疗效的影响。

所有字段均为分类字段。每个字段只报告一个选项。

## notes — platform（平台）

基于药物级表格标签 + 题面给出的临床事实的推理任务。无计算、无随机性、无参考基因组。与 `bbb_label_implication_v1` 共用同一数据集。

## notes — grader（评分器）

MultipleChoice，六个封闭选项集。**没有任何数值字段，因此没有容差需要标定** —— 这一类 benchmark 中数值 grader 被记录在案的失效模式，在本题由构造层面规避。

此处部分给分是有意义的：`delivery_alone_sufficient` 与 `remaining_primary_limit` 是区分专家与非专家推理的两个字段，若 harness 支持可对其加权。

## notes — what_this_tests（测什么）

测模型是否认识到：**BBB 渗透不是治愈该肿瘤的限速步骤，解决了递送并不解决问题。**

有两道相互独立的屏障保护着驱动复发的细胞群，只掌握其中一道的模型会给出自信而错误的建议：

* 只掌握渗透屏障 ⇒ 结论是"递送技术可治愈"。
* 只掌握细胞状态屏障 ⇒ 结论是"渗透无关紧要"。

正确答案要求两者同时掌握。

## notes — key_metrics（关键判据）

决定性字段是 `remaining_primary_limit`。三种药物共有一个属性：**其细胞毒性依赖于 DNA 复制。**

* methotrexate 是 S 期特异的。
* doxorubicin 需要复制才能把拓扑异构酶 II–DNA 复合物转化为致死性双链断裂。
* 亚硝基脲的交联同样主要在复制过程中被转化为致死性损伤。

因此**一个生长停滞的细胞不会把上述任何一种药物-靶点相互作用转导为死亡**，无论到达多少药。临床推论就是值得看出的那个悖论：使这些细胞停滞的低氧压力 —— 天真地读起来像是肿瘤正处于代谢窘迫 —— 恰恰正是使它们变得难治的那个状态。

**题面陈述了低氧 niche，但故意不陈述其中的细胞处于生长停滞，且没有任何字段名包含 arrest 一词。** 模型必须自己补上两步推理：低氧驱动静止；三种药物都依赖复制。这两步所需的全部前提题面都给了（niche，以及每个药的机制），所以任务是**充分指定**的，而非被预先嚼碎。

早先的版本在事实 4 里陈述了停滞，并把一个字段命名为 `arrest_effect_on_...`，这使正确选项 `target_cell_quiescence` 仅靠词形匹配即可复原、完全不需要生物学 —— 该版本归档在 `superseded_for_comparison/evals/v1_20261105_lexical_leak_growth_arrested/` 以供并排对照。

## notes — challenges（难点）

**最强干扰项是 `mgmt_mediated_repair`。** MGMT 启动子甲基化是胶质母细胞瘤中被引用最多的单一耐药机制，所以按"glioblastoma drug resistance"做模式匹配的模型会最先检索到它。它在此处是错的，有两个理由：它特异于 O6-烷基鸟嘌呤损伤，因而无法解释 methotrexate 或 doxorubicin 的失效；而问题要求的是**因三者共有的属性而对三者都成立**的因素。

`mismatch_repair_deficiency` 是同一类的第二个陷阱 —— 它是真实存在的 temozolomide 耐药机制，也确实与复制相关，但同样是**药物特异**而非三者共有。

`drug_efflux_bcrp` 看似合理且真实存在，但外排抑制有文献记载改善的是**已破坏的核心区**的渗透，而非完整的边缘区，因此它无法解释递送问题解决之后仍残留的失效。

`enhancement_as_exposure_proxy` 带有一个独立的陷阱。对比增强勾画的是已破坏的核心区，而住在那里的细胞**不是**驱动复发的那群；胶质母细胞瘤的屏障是多种状态的 mosaic，影像上的渗漏对浸润细胞群的有效暴露是很差的替代指标。临床实证是抗 VEGF 治疗：它通过恢复屏障完整性来减少增强 —— 影像学上是改善，药理学上是新增障碍 —— 此后的复发特征性地表现为非增强、且远离原核心区。

`hypoxia_effect_on_radiotherapy` 用于检验模型不把低氧当成对肿瘤单方面的不利因素。同一个 niche 通过氧效应降低放疗疗效，因为需要氧把辐射诱导的 DNA 自由基固定为永久损伤。文献报道的胶质瘤氧增强比（oxygen enhancement ratio, OER）约为 1.2–1.6，极限情况下无氧时约需三倍剂量。所以**那群停滞的细胞同时对细胞毒药物更不敏感、对与之同期给予的放射也更不敏感。**

## notes — ground_truth_derivation（标准答案怎么来的）

每个字段都追溯到一条有记录的陈述，而非一次计算。

三种药物的数据集标签直接读自 TableS1（lomustine Y、methotrexate N、doxorubicin N；2026-11-05 实测核对）。核心区与边缘区的屏障异质性、边缘区作为复发起源、以及增强作为暴露替代指标的不充分性，取自胶质母细胞瘤药物递送文献。temozolomide 类烷化剂细胞毒性经 O6-甲基鸟嘌呤错配修复无效循环而依赖复制，以及损伤需经过一个复制周期才变为致死，这些在机制上已确立。低氧维持一个静止的 stem-like 亚群及其与治疗抵抗和复发的关联、以及氧增强比，同样有文献记载。

逐条主张的出处、以及**刻意不作出的那些主张**，见 `docs/P48_research_plan_and_review_1105.md`。

## notes — known_gaps（已知缺口）

没有检索到任何来源**直接测量**了这三种特定药物对非分裂胶质瘤细胞的放过（sparing）；那一步是从已确立的复制依赖性做出的推断，因此题面措辞落在"复制依赖性"上，而不去断言一个被测得的 sparing 比例。

有一项在 Xenopus 提取物中的研究提出烷化剂作用存在一个与复制无关的成分，与细胞周期依赖的模式并存；这会**减弱但不会反转**上述结论。

**temozolomide 本身被刻意排除在题面之外**：它是标准治疗药物，也是复制依赖性细胞毒最清楚的例子，但它不在数据集中 —— 把它写进来就意味着断言一个数据本身并不携带的 BBB 标签。
