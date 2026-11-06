# P48 — LatchBio「Question Creation」实际题面（2026-11-05 由截图转录）

> **来源**：用户 2026-11-05 提供的 4 张测试界面截图。**这是对方的真实要求**，
> 此前项目里所有设计依据的是招聘邮件 `Email.txt:19` 那一句概述，**与本文有实质差异**。
> 本文为转录，非逐字扫描；关键约束已逐条核对。若与界面原文冲突，以界面为准。

## Background（对方的立意）

> At Latch, we build problems, to understand and improve how large language models reason
> about **oncology**. A good problem gives an LLM **the evidence it needs** and asks it to make
> a defensible scientific analysis or decision. It should be challenging **because of the
> scientific judgment involved, not because important information is missing or the
> instructions are unclear.**

⇒ **"难"必须来自科学判断，不能来自藏信息。** 这与本项目 §3.7.2 自己总结的规格原则方向一致，
但对方把界线划得更靠"给足证据"一侧。

## Your task

选**一篇**肿瘤学领域的原始研究论文（primary research paper in oncology，可以是基础癌生物学、
转化研究或临床研究），**来自你有直接研究或应用经验的方向**。论文应围绕一个肿瘤学问题，
例如 **biomarker identification、patient stratification、therapeutic regimen efficacy**。

### ⛔ 三条硬性限制

| 限制 | 原文 |
| :--- | :--- |
| **不能是自己的论文** | `The paper must not be your own, and you should not be one of its authors.` |
| 只能用正式发表版 | `Use the final published version` |
| 只能用公开信息 | `and only publicly available information` |

论文**不必新、也不必有名**。

## 1. Why this paper?

给出**引文 + 链接**，后接 **3–4 条简洁 bullet**：

- 这篇论文为什么留在你心里
- 其中某一个想法 / 结果 / 决定如何改变了你思考肿瘤学的方式
- 它的影响在你自己的科研工作中体现在哪里
- 一条你至今认为重要的 insight、局限或未决问题

> `Do not summarize the paper or restate its abstract. We value specific, insightful points
> over a long description.`

## 2. Revisit one result

选论文中**一张图、一张表或一个核心结果**，描述：

- 作者想回答的问题
- 他们得出的结论
- **如果你是作者，你会改变分析或解释中的什么**
- 为什么你提出的改变重要
- 它会**加强、削弱还是实质性改变**原结论

> 若公开数据不足以做该分析，就**说明你会要求做哪一项确切的额外分析或对照**，
> 以及**什么样的结果能区分相互竞争的解释**。

## 3. Create, test, and revise two LLM problems

两道题，考科学推理，在你有直接经验的肿瘤学方向。

| 约束 | 要求 |
| :--- | :--- |
| 与第 2 部分的关系 | **至少一道必须直接从你在第 2 部分讨论的那个结果生长出来**；第二道可用同一篇论文的另一个结果 |
| 题目性质 | 必须反映**真实的研发判断**：选对照、选比较、选 metric、选阈值、assay 解读、推理分析、生物学结论 |
| 禁止 | **不得是** paper-summary、literature-lookup、或能从标题/摘要生成的 fact-recall |
| 定量 | **至少一道必须要求定量或基于数据的分析**；另一道可聚焦生物学解释或决策 |
| ⛔ 自足性 | `Both must be answerable using only the materials you provide; the LLM should not need to browse the paper or recognize it.` |

### 每道题必须包含的七个部分

| # | 部分 | 要求 |
| :---: | :--- | :--- |
| 1 | **Task prompt** | 你要给 LLM 的**确切 prompt** |
| 2 | **Inputs** | 提供给 LLM 的数据/图/表/摘录。**注明来源并说明你做过的任何处理**。不得含机密或专有信息 |
| 3 | **Expected answer** | 给出答案，并**展示它如何由所提供的证据推导出来** |
| 4 | **Grading rule** | 定义结构化输出与**确定性**评分规则。须写明：哪些字段/数值被评分；可接受答案与数值容差；**科学上同样站得住、也应判通过的替代答案**；什么条件导致判不通过 |
| 5 | **Scientific skill** | 这道评测测的是什么、为什么在实践中重要、**你的领域经验如何塑造了它** |
| 6 | **Plausible wrong approach** | 至少一条看起来合理、却会得出不同答案的捷径或错误。评测必须能把它与站得住的分析区分开，**且不能依赖信息缺失、措辞含糊或格式把戏** |
| 7 | **LLM testing and revision** | 见下 |

> 评分脚本可以附，但**一条精确到别人无需再做科学判断即可实现的规则或伪代码就够**。
> **评的是最终的科学答案，不是文风。**

### 第 7 部分：LLM 实测与修订（**硬性，不可省**）

**至少用一个常见 LLM 测两次**（免费账号即可）。必须包含：

- 用的**模型/产品名与测试日期**
- **完整响应**，放进附录
- 对每条响应**应用你的评分规则的结果**
- 模型哪里做得好、推理在哪里失败
- 测试后你对 prompt / inputs / expected answer / grader 做的**任何改动**

> 若你认为无需修订，**要解释为什么**。
> `We are interested in how you use LLM behavior to improve the evaluation, not simply
> whether the model passes or fails.`

## Example Problem（对方给的范例，来自 benchmarks.bio）

范例是 GBM 单细胞 RNA-seq 的一道题（给 raw counts + cell metadata + primary-target 注释，
问哪个靶点在脑/肿瘤 compartment 相对免疫细胞最一致地富集，报告 specificity score，再选解释）。

**⭐ 关键是它的输出格式**，我们之前反推的 `eval_answer.json` 是错的：

```
Return EXACTLY:

<EVAL_ANSWER>
{"top_compartment_specific_target": ["<gene>"], "compartment_specificity_score": <float>, "answer": "<letter>"}
</EVAL_ANSWER>

Report the score rounded to two decimals, where <letter> is A, B, C, or D.
Provide a one-line rationale outside the answer tags.
```

⇒ **答案包在 `<EVAL_ANSWER>` 标签里**，标签外给一行 rationale。范例还显示：
多选题用 A/B/C/D 字母、数值明确要求保留两位小数。

## What to submit

**一份 Markdown 或 PDF 文档**，外加任何支撑数据与评分代码。文档分五节：

| # | 节 |
| :---: | :--- |
| 1 | Paper choice |
| 2 | Result critique |
| 3 | Problem 1 |
| 4 | Problem 2 |
| 5 | LLM outputs appendix |

> 可以用 LLM 完成本练习，但要**简要说明用了哪些工具、怎么用的**。
> **科学选择、预期答案与评分决定必须是你自己的。**

## What we are looking for

- 肿瘤学 / 癌症研究的知识深度
- 对已发表证据的深思熟虑的分析
- **扎根于真实科研工作**的问题，而非泛泛的复杂度
- 公平、自足、答案站得住的评测
- 清晰且确定性的评分
- **能测 LLM、诊断其行为、并据此改进评测的证据**
- 简洁、直接的科学表达

## ⚠ 与本项目既有设计的差距

见 `P48_research_plan_and_review_1105.md` §5。
