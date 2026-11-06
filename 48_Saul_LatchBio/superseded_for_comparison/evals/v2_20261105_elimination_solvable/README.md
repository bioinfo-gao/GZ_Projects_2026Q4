# v2 归档 — `glioma_margin_rate_limit_v1` 的「排除法可解」版

归档日期：2026-11-05
归档原因：**题目不需要 glioma 知识即可全对。决定性字段可用纯排除法解出，且题面自带提示句指向那条排除路径。**

替代版本：`evals/glioma_margin_rate_limit_v3.json`（id 改为 `_v3`，陷阱已更换，不是措辞修订）。

## 本目录内容

| 文件 | 说明 |
| :--- | :--- |
| `glioma_margin_rate_limit_v1.json` | **v2 原件，逐字节原样**（`cmp` 验证与被替换前的现行件相同） |
| `rendered_v2_for_reading.md` | v2 的英文可读渲染 |
| `rendered_v2_zh.md` | v2 的中文审阅译文 |
| `README.md` | 本文 |

**这是原始产物本身，不是事后复现。** v2 在被 v3 取代前先 `cp -p` 到此处，时间戳保留。
⚠ 文件名仍是 `_v1`：v1→v2 是同一文件原地修订（只改了字面泄题），v2 才拆出独立 id。
上一层 `v1_20261105_lexical_leak_growth_arrested/` 里的同名文件是**更早的 v1**，两者不同。

## 缺陷是什么

### 缺陷 1（主）：决定性字段可纯排除法解出

`remaining_primary_limit` 的五个选项里，四个干扰项**都是药物特异**的，只有正确选项是药物无关的：

| 选项 | 不懂 glioma 也能排掉的理由 |
| :--- | :--- |
| `insufficient_dose` | 题面前提已声明"与核心区同浓度" |
| `drug_efflux_bcrp` | 外排属于递送，前提已声明递送解决 |
| `mgmt_mediated_repair` | 通识即知只管烷化剂，三药不共有 |
| `mismatch_repair_deficiency` | 同上，TMZ 特异 |
| `target_cell_quiescence` | **剩下的唯一一个** |

而题面里我亲手写了一句 `Choose the factor that applies to all three agents by virtue of a property they share` ——
**这句话等于把排除法的判据直接交给模型。** 它本意是"充分指定"，实际是把唯一的药物无关选项点了出来。

### 缺陷 2：六个字段里四个是白送或通识

| 字段 | 答对需要什么 | 判定 |
| :--- | :--- | :---: |
| `agent_reaching_margin_unaided` | 读题面给的标签 | 白送 |
| `enhancement_as_exposure_proxy` | 事实 1–3 已把答案摆好 | 白送 |
| `delivery_alone_sufficient` | benchmark 里这种问法基本不是 yes | 强暗示 |
| `remaining_primary_limit` | 见缺陷 1 | 排除法可解 |
| `hypoxia_effect_on_cytotoxic_sensitivity` | 低氧→化疗抗性，教科书 | 通识 |
| `hypoxia_effect_on_radiotherapy` | 氧效应 / OER，教科书 | 通识 |

另外最后两个字段近乎同一问题的两问（都是"低氧对治疗不利"），不提供独立鉴别力。

### 缺陷 3（根因）：陷阱方向与文献主流叙事一致

低氧 → 停止增殖 → 复制依赖药杀不死，这条链在 glioma 综述里是**标准叙事**。
**凡是正确答案与主流叙事同向的陷阱都不可能难** —— 模型检索到主流叙事就对了。
真正难的陷阱必须满足：**主流叙事给出的答案是错的。**

### 缺陷 4（用户 2026-11-05 指出）：把 BBB 当成了真实障碍，而它应该是诱骗项

> 「BBB 在脑瘤中不是特别重要，因为 glioma 绝大多数伴随 BBB 破坏，可以作为诱骗项。」

v2 的整个前提是"穿不过 BBB 所以到不了"。但增强区 BBB 早已自我破坏，药物的完整-BBB 通透性
在那里并不决定暴露；而在屏障完整的边缘区，暴露又不决定杀伤（细胞不增殖）。
⇒ **`BBB +/-` 标签在这个病里几乎不决定结局，它是数据集里最显眼、却不是起作用的那个变量。**
v2 没有利用这一点，反而把它当作真实限速步骤，等于放弃了最好的陷阱。

## v3 怎么改的

| 改动 | 针对 |
| :--- | :--- |
| 删掉 `...by virtue of a property they share` 提示句 | 缺陷 1 |
| 决定性字段五个选项**全部改为药物无关**（屏障通透性 / 增殖状态 / 血管外间质转运 / DNA 修复能力 / 免疫微环境） | 缺陷 1 |
| 删掉两个白送字段，合并两个重复的 hypoxia 字段，6 字段 → 5 字段 | 缺陷 2 |
| 新增两个**方向反转**字段：持续抗 VEGF 对边缘区暴露、对内在敏感性的影响，正确答案与 vascular normalization 主流框架相反 | 缺陷 3 |
| 新增 `agent_selection_by_bbb_label`：承认标签上 lomustine 占优、同时否认它预测疗效 —— 把 BBB 做成诱骗项 | 缺陷 4 |
| 新增 `relapse_pattern_after_antivegf` 作为**自洽性检验**：方向答反的模型在此处会给出自己的机制无法产生的复发模式 | 缺陷 3 |

## 怎么并排对照

```bash
cd /home/gao/projects_2026Q4/48_Saul_LatchBio
A=superseded_for_comparison/evals/v2_20261105_elimination_solvable

# 1) 题面逐行对照（v2 vs v3）
diff <(python3 -I -c "import json;print(json.load(open('$A/glioma_margin_rate_limit_v1.json'))['task'])") \
     <(python3 -I -c "import json;print(json.load(open('evals/glioma_margin_rate_limit_v3.json'))['task'])")

# 2) 两版答案字段并列
python3 -I -c "
import json
for p in ['$A/glioma_margin_rate_limit_v1.json','evals/glioma_margin_rate_limit_v3.json']:
    d=json.load(open(p)); print('==',d['id'])
    for k,v in d['grader']['config']['ground_truth'].items(): print('  %-48s %s'%(k,v))
"

# 3) 读旧版可读渲染（中/英各一份；新版渲染在 tmp/）
cat "$A/rendered_v2_zh.md"
```

## 复查判据

对照时该确认的是：**v3 的决定性字段是否真的无排除法路径。**
机械检查 = 列出 `primary_determinant_of_margin_response` 的五个选项，逐个问
「这个选项是否因为只适用于三药中的某一个而能被排掉」，五个都应回答"否"。

```bash
cd /home/gao/projects_2026Q4/48_Saul_LatchBio
python3 -I -c "
import json,re
t=json.load(open('evals/glioma_margin_rate_limit_v3.json'))['task']
print(re.search(r'primary_determinant_of_margin_response. \(string\) — one of (.+?)\. The',t).group(1))
"
```

2026-11-05 实测：五个选项分别是屏障通透性、靶细胞增殖状态、血管外间质转运、DNA 修复能力、
免疫微环境 —— 全部是肿瘤或递送侧的属性，没有一个是"某一个药的属性"，排除法失效。

## 尚未解决

**v3 同样未在任何模型上实测过。** v2 的缺陷是我自己审计出来的，不是实测出来的；
v3 是否真的更难，仍需实测证实。两个方向反转字段的风险在于：若某模型恰好读过 GBM 里
抗 VEGF 屏障恢复那一支文献，它会直接答对，此时难度来自知识覆盖而非推理，
仍需考虑进一步加深。

归档内容一律不进交付目录。
