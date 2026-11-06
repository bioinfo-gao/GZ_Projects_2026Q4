# glioma_margin_rate_limit_v3（中文审阅版）

> ⚠ **本文是译文，不是真身。** 内容真身在 `evals/glioma_margin_rate_limit_v3.json`（英文）。
> 要改内容请改那个 JSON，然后重译本文。交付给 LatchBio 的是英文 JSON，本文不交付。
> 翻译日期：2026-11-05 ｜ v2 归档在 `superseded_for_comparison/evals/v2_20261105_elimination_solvable/`

**数据：** https://github.com/bioinformatics-gao/CASE-BBB-prediction-Data

## Grader — MultipleChoice

| 字段 | 标准答案 |
| :--- | :---: |
| `primary_determinant_of_margin_response` | `proliferative_state_of_the_target_cells` |
| `agent_selection_by_bbb_label` | `lomustine_preferred_but_this_does_not_predict_response` |
| `antivegf_effect_on_margin_exposure` | `decreases_via_barrier_restoration` |
| `antivegf_effect_on_margin_cytotoxic_sensitivity` | `decreases` |
| `relapse_pattern_after_antivegf` | `non_enhancing_and_distant` |

全部精确匹配，无容差可标定。

## 题面（交给模型的原文，此处为中译）

一位患者患复发性胶质母细胞瘤（glioblastoma）。有三种细胞毒药物可用，其血脑屏障（blood-brain barrier, BBB）通透性标签记录在随附数据集中（`1) Supplementary-data/TableS1-Drug-for-Train-and-Validation.xlsx`，列 `BBB +/-`）：

* **lomustine** —— 数据集标签 `Y`；亚硝基脲类，烷基化并交联 DNA
* **methotrexate** —— 数据集标签 `N`；抗叶酸药，阻断胸苷酸合成，S 期特异
* **doxorubicin** —— 数据集标签 `N`；拓扑异构酶 II 毒物，在复制过程中把酶-DNA 复合物转化为双链断裂

该肿瘤具有以下事实：

1. 对比增强的肿瘤主体 BBB 已被破坏，血管通透性升高。
2. 非增强的浸润性边缘区 BBB 性质基本完整。
3. 复发源自边缘区的浸润细胞，而非主体。
4. 边缘区含有一个低氧（hypoxic）niche。

此外正在考虑加用贝伐珠单抗（bevacizumab，抗 VEGF）。**请按持续抗 VEGF 治疗作答，即超出任何短暂的血管正常化（vascular-normalisation）窗口之后。**

把答案以 JSON 写入 `eval_answer.json`，字段恰好为：

* `primary_determinant_of_margin_response`（字符串）—— `bbb_permeability_of_the_agent`、`proliferative_state_of_the_target_cells`、`interstitial_transport_beyond_the_vessel_wall`、`agent_specific_dna_repair_capacity`、`immunosuppressive_microenvironment` 之一。主要决定这些药物能否减少驱动复发的细胞群的那个因素。
* `agent_selection_by_bbb_label`（字符串）—— `lomustine_preferred_and_this_predicts_response`、`lomustine_preferred_but_this_does_not_predict_response`、`no_label_based_preference_exists`、`methotrexate_preferred_and_this_predicts_response` 之一。数据集的 `BBB +/-` 标签对这位患者在三药之间的选择有何影响。
* `antivegf_effect_on_margin_exposure`（字符串）—— `increases_via_vascular_normalisation`、`increases_via_reduced_interstitial_pressure`、`decreases_via_barrier_restoration`、`unchanged_because_margin_barrier_already_intact` 之一。持续抗 VEGF 治疗对浸润性边缘区细胞毒药物暴露的影响。
* `antivegf_effect_on_margin_cytotoxic_sensitivity`（字符串）—— `increases`、`decreases`、`no_effect` 之一。持续抗 VEGF 治疗对边缘区细胞群**对这三种药物的内在敏感性**的影响；**不计**递送药量本身的变化。
* `relapse_pattern_after_antivegf`（字符串）—— `enhancing_and_local`、`non_enhancing_and_distant`、`unchanged_from_baseline_pattern` 之一。持续抗 VEGF 治疗下复发的特征性影像学模式。

所有字段均为分类字段。每个字段只报告一个选项。

## notes — platform（平台）

基于药物级表格标签 + 题面给出的临床事实的推理任务。无计算、无随机性、无参考基因组。

与 `bbb_label_implication_v1` 共用同一数据集，且**是刻意设计的反向配对**：那道题要求把 BBB 标签当作数据认真对待，这道题要求认识到这些标签并不支配临床结局。

## notes — grader（评分器）

MultipleChoice，五个封闭选项集。无数值字段，因此无容差需要标定 —— 这一类 benchmark 中数值 grader 被记录在案的失效模式，在本题由构造层面规避。

五个字段**彼此独立鉴别**，没有一个是另一个的改写。其中两个（`antivegf_effect_on_margin_exposure` 与 `antivegf_effect_on_margin_cytotoxic_sensitivity`）的正确答案**与主流发表框架相反**，因此无法通过检索共识语句得到。

## notes — what_this_tests（测什么）

**测模型能否识别出：数据里最显眼的变量不是起作用的那个。**

这是一个 BBB 通透性数据集，题面又以 BBB 标签开头，所以牵引力全部指向"把穿过 BBB 当作决策变量"。而在胶质母细胞瘤里它基本不是：

* 对比增强的主体**已经自行破坏了屏障**，所以某药在完整 BBB 下的通透性并不决定它在那里的暴露；
* 而在屏障完整的边缘区，暴露**仍然不决定杀伤**，因为靶细胞并不处于增殖状态。

只按标签推理的模型会给出一个自信的、有数据支撑的、错误的建议 —— 它会偏向 lomustine，并预期这个偏向会转化为疗效。

后半部分测的是：模型持有的是**胶质母细胞瘤里抗 VEGF 效应的方向**，还是血管正常化文献所断言的方向。这两支文献彼此矛盾，而对 benchmark 有意义的问题是模型把哪一支套用到这个病上。

## notes — key_metrics（关键判据）

**决定性的是 `primary_determinant_of_margin_response` 与 `agent_selection_by_bbb_label` 这一对。**

模型完全可以把表读对 —— lomustine 是唯一的 `Y` —— 却仍然答错，因为问的是**这个读法能不能预测疗效**。正确答案**承认标签上的偏向、同时否认它的预测力**，这是唯一一个要求同时掌握两半图景的组合。

`no_label_based_preference_exists` 在数据上就是错的：标签确实有差别，lomustine 确实是能到达完整屏障的那一个。

否认的机制在于三种药物**都依赖复制**：methotrexate 是 S 期特异；doxorubicin 需要复制才能把拓扑异构酶 II–DNA 复合物转化为致死性双链断裂；亚硝基脲的交联同样主要在复制过程中被转化为致死性损伤。因此一个不在增殖的细胞不会把上述任何药物-靶点相互作用转导为死亡，无论到达多少药。

**题面陈述了低氧 niche，但不陈述其中细胞不在增殖，且没有任何字段名包含 arrest 或 quiescence。** 两步推理都留给模型：低氧抑制增殖；三药都依赖复制。所需前提题面全部给了。

## notes — challenges（难点）

四个陷阱彼此独立。

**第一，BBB 诱骗项。** 随附数据集是 BBB 通透性数据集，每个药都以其 BBB 标签引入，开头两条事实都在讲屏障状态。答 `bbb_permeability_of_the_agent` 的模型**不是在草率地模式匹配 —— 它是在依据所给数据最显著的特征进行推理**，这正是这个诱骗项强的原因。

该字段的选项集**不留排除法路径**：五个选项全是肿瘤侧或递送侧的属性，而不是某一个药的属性，所以解掉本题早先版本的那条"共有属性"捷径在此失效。`interstitial_transport_beyond_the_vessel_wall` 与 `immunosuppressive_microenvironment` 都真实且被大量引用；`agent_specific_dna_repair_capacity` 是 MGMT 那个答案的一般化形式。

**第二，抗 VEGF 对暴露影响的方向。** 四个选项里有三个用的是血管正常化框架自己的词汇，而 `increases_via_reduced_interstitial_pressure` 是一个真实且被广泛发表的机制。在胶质母细胞瘤里，持续抗 VEGF 治疗**恢复**屏障完整性，从而**减少**药物进入非增强边缘区 —— 方向相反。增强减少正是这种屏障恢复的影像学特征，这也是它读起来像好转的原因。

**第三，抗 VEGF 对内在敏感性影响的方向。** 把抗 VEGF 与细胞毒药物、放疗联用的已发表理由就是"正常化改善氧合、因而提高敏感性"。而在持续治疗下，血管修剪反而**加深**边缘区低氧，进一步抑制增殖，于是**降低**对三种复制依赖药物的敏感性。该字段明确要求"不计递送药量变化"，所以模型不能靠复用上一个字段的答案作答，必须真的把低氧→增殖→复制依赖这条链走完。

**第四，`relapse_pattern_after_antivegf` 是自洽性检验。** 若模型把前两个方向字段按正常化框架作答、却在这里答 `non_enhancing_and_distant`，它陈述的复发模式是它自己给出的机制无法产生的。因此**联合评分**可以区分"持有一套机制"的模型与"检索到若干事实"的模型。

## notes — ground_truth_derivation（标准答案怎么来的）

每个字段都追溯到一条有记录的陈述，而非一次计算。

三药的数据集标签直接读自 TableS1（lomustine Y、methotrexate N、doxorubicin N；2026-11-05 实测核对）。增强主体内屏障破坏 vs 非增强浸润边缘区屏障相对完整、以及边缘区作为复发起源，取自胶质母细胞瘤药物递送文献。烷化剂细胞毒性经 O6-甲基鸟嘌呤错配修复无效循环而依赖复制、以及损伤需经过一个复制周期才变致死，在机制上已确立。低氧维持一个不增殖的 stem-like 亚群及其与治疗抵抗和复发的关联，同样有记载。抗 VEGF 恢复屏障完整性伴增强减少、随后特征性的非增强且远处复发模式、以及持续治疗下血管修剪伴肿瘤低氧加深，是已报道的胶质母细胞瘤观察。

逐条主张的出处与**刻意不作出的主张**，见 `docs/P48_research_plan_and_review_1105.md`。

## notes — known_gaps（已知缺口）

三个缺口，都不改变答案，但都界定了主张的边界。

**一、血管正常化窗口确实存在。** 有记载的短暂正常化窗口期内，氧合与递送确实改善。因此题面在前言和每个抗 VEGF 字段的描述里都**明确限定为该窗口之后的持续治疗**。若不作此限定，`increases_via_reduced_interstitial_pressure` 就是站得住的，该题会变成**有歧义**而不是难。

**二、持续抗 VEGF 下低氧加深更多依赖临床前模型**，而非在浸润边缘区组织里的直接人体测量；临床推断来自复发表型，不是来自边缘区的氧电极。答案 `decreases` 是站得住的方向，但它是推断，不是测量。

**三、没有任何检索到的来源直接测量了这三种特定药物对不增殖胶质瘤细胞的放过（sparing）。** 那一步是从已确立的复制依赖性做出的推断，因此推理落在"复制依赖性"而非一个被断言的 sparing 比例上。有一项在 Xenopus 提取物中的研究提出烷化剂作用存在一个与复制无关的成分、与细胞周期依赖模式并存，这会减弱但不会反转结论。

**temozolomide 被刻意排除**：它是标准治疗药物，也是复制依赖性细胞毒最清楚的例子，但它不在数据集中 —— 写进来就意味着断言一个数据本身并不携带的 BBB 标签。

## notes — version_history（版本沿革）

v1 与 v2 归档在 `superseded_for_comparison/evals/` 下。

* **v1** 在事实 4 里陈述了边缘区细胞生长停滞，并把一个字段命名为 `arrest_effect_on_...`，使正确选项仅靠词形匹配即可复原、完全不需要生物学。
* **v2** 去掉了那个字面泄题，但**仍可在不懂 glioma 的情况下解出**：六个字段里四个是题面白送或教科书通识；决定性字段可用排除法解出，因为四个干扰项都是药物特异、而正确选项是唯一药物无关的那个 —— 而题面里那句"选择因三药共有属性而成立的因素"**主动给这条排除路径立了路标**。
* **v3** 删掉该提示句，把该字段的所有选项都改为药物无关，删掉两个可由读题面得到的字段，合并两个重复的 hypoxia 字段，并换成两个正确答案与主流发表框架相反的方向反转字段。
