# v1 归档 — `glioma_margin_rate_limit_v1` 的字面泄题版

归档日期：2026-11-05
归档原因：**题面存在字面泄题，正确选项可仅靠词形匹配选出，不需要任何生物学推理。**

## 本目录内容

| 文件 | 说明 |
| :--- | :--- |
| `glioma_margin_rate_limit_v1.json` | **v1 原件，逐字节原样**（`cmp` 验证与当时的现行件相同） |
| `rendered_v1_for_reading.md` | v1 的人类可读渲染（由 `scripts/2_render_eval.py --md` 产生） |
| `README.md` | 本文 |

**这是原始产物本身，不是事后复现。** v1 在被覆盖前先 `cp -p` 到此处，时间戳保留。

## 缺陷是什么

v1 的题面把答案写在了premise 里：

| 位置 | v1 原文 |
| :--- | :--- |
| 临床事实 #4 | `The margin contains a hypoxic niche in which a substantial fraction of tumour cells are **growth-arrested**.` |
| 答案字段名 | `**arrest**_effect_on_cytotoxic_sensitivity` |
| 正确选项 | `target_cell_**quiescence**` |

`growth-arrested` 与 `quiescence` 近义，且字段名里也有 `arrest`。
⇒ **模型只要在题面里找与选项同义的词，就能选中 `target_cell_quiescence`，
完全不需要知道"低氧→静止"或"三药都依赖 DNA 复制"。**
这使本题退化为词形匹配题，测不出设计意图中的科学深度。

发现经过：把 JSON 渲染成可读形式后，题面与选项并排出现，线索才变得显眼。
raw JSON 里这两处相距很远且被 `\n` 转义包裹，不易察觉。

## 怎么改的（用户 2026-11-05 选定方案 A）

| 处 | v1 | v2（现行） |
| :--- | :--- | :--- |
| 事实 #4 | `a hypoxic niche in which a substantial fraction of tumour cells are growth-arrested.` | `a hypoxic niche.` |
| 字段名 | `arrest_effect_on_cytotoxic_sensitivity` | `hypoxia_effect_on_cytotoxic_sensitivity` |
| 该字段描述 | `The effect of hypoxia-induced growth arrest on...` | `The effect of the hypoxic state at the margin on...` |

**改后 `arrest` 一词在全题面消失**，`quiescence` 只作为 5 个选项之一出现
（必须保留，否则无法作答）。

**为什么这样改仍然是"充分指定"而非"故意藏信息"**：作答所需的全部前提都仍在题面 ——
低氧 niche（事实 #4）、以及每个药的机制（S 期特异 / 需复制）。
删掉的只是**本该由作答者自己完成的那一步推理**。
LatchBio `adding_evals.md` 要求任务描述指明方法、约束与输出格式，这三者未受影响。

另一个被否的方案 B 是"保留事实 #4、只改选项名"，否决理由：
"静止"这个概念仍已出现在题面，模型仍可据此缩小选项范围。

## 怎么并排对照

```bash
cd /home/gao/projects_2026Q4/48_Saul_LatchBio
A=superseded_for_comparison/evals/v1_20261105_lexical_leak_growth_arrested

# 1) 结构化 diff：只看三处改动
diff <(python3 -I -m json.tool "$A/glioma_margin_rate_limit_v1.json") \
     <(python3 -I -m json.tool evals/glioma_margin_rate_limit_v1.json)

# 2) 题面逐行对照（推荐，最直观）
diff <(python3 -I -c "import json;print(json.load(open('$A/glioma_margin_rate_limit_v1.json'))['task'])") \
     <(python3 -I -c "import json;print(json.load(open('evals/glioma_margin_rate_limit_v1.json'))['task'])")

# 3) 读旧版可读渲染（新版渲染在 tmp/）
cat "$A/rendered_v1_for_reading.md"
```

## 复查判据

对照时该确认的是：**v2 题面里是否还存在任何与 `target_cell_quiescence` 近义的词**。
机械检查方式：

```bash
cd /home/gao/projects_2026Q4/48_Saul_LatchBio
python3 -I -c "
import json,re
t=json.load(open('evals/glioma_margin_rate_limit_v1.json'))['task']
for w in ['arrest','quiescen','dormant','non-dividing','non-proliferat','senescen','G0','cell.cycle']:
    print('%-16s'%w,[m.start() for m in re.finditer(w,t,re.I)])
"
```

2026-11-05 实测：除 `quiescen` 在选项列表中出现 1 次（必须保留）外，其余全部为空。

## 尚未解决

**v2 仍未在任何模型上实测过。** 修掉字面泄题只是排除了一条平凡捷径，
不等于强模型答不出来 —— 这条推理链在 glioma 文献里并不罕见。
实测结果若显示强模型仍直接答对，则需更换陷阱，而不是继续打磨措辞。

归档内容一律不进交付目录。
