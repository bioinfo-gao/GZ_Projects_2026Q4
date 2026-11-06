# v3 归档 — 「前提陈述了该考的生物学」版

归档日期：2026-11-05
归档原因：**题面用四条"临床事实"把答案所需的疾病生物学直接给了模型。其中第 1 条恰恰就是让 BBB 标签成为诱骗项的那条知识。**

替代版本：`evals/glioma_margin_rate_limit_v4.json`

## 本目录内容

| 文件 | 说明 |
| :--- | :--- |
| `glioma_margin_rate_limit_v3.json` | **v3 原件，逐字节原样**（`cmp` 验证） |
| `rendered_v3_for_reading.md` | v3 英文可读渲染 |
| `rendered_v3_zh.md` | v3 中文审阅译文 |
| `README.md` | 本文 |

## 缺陷是什么（用户 2026-11-05 提问指出）

> 「为什么要告诉 AI 这些事实？ 该肿瘤具有以下事实：对比增强的肿瘤主体 BBB 已被破坏，血管通透性升高。」

v3 的题面有一个四条事实的清单。**这四条全是胶质母细胞瘤的一般生物学，不是这个患者的情况**：

| v3 的"事实" | 性质 | 本应 |
| :--- | :---: | :---: |
| 1. 增强主体 BBB 已破坏、通透性升高 | GBM 通识 | **删** |
| 2. 非增强浸润边缘屏障基本完整 | GBM 通识 | **删** |
| 3. 复发源自边缘浸润细胞而非主体 | GBM 通识 | **删** |
| 4. 边缘区含低氧 niche | GBM 通识 | **删** |

**第 1 条是致命的那一条。** v3 的整个立意是"BBB 标签是诱骗项"，而诱骗项之所以有效，
**全靠模型自己不知道标签不管用**。写上第 1 条等于直接告诉它标签不管用 ——
陷阱被自己的题面拆除了。

v3 还逐药写了机制（`S-phase specific`、`converts ... during replication`、`alkylates and crosslinks DNA`），
那是标准药理学，等于供出了**复制依赖性**这一半推理。

## 根因：把"问题充分指定"误当成"前提喂到嘴里"

LatchBio 的 `adding_evals.md` 要求任务描述指明方法、约束与输出格式。
我据此把生物学也一并写了进去。**正确的界线是：**

> **患者层面、因病例而异、无法先验推知的 ⇒ 必须写明。**
> **疾病层面、对所有 GBM 都成立的 ⇒ 一律不写，那正是考题。**

⇒ 该原则已写入 v4 的 `notes.specification_principle`，成为本题的显式规格。

## ⚠ 这是同一个错误的第三次

| 版本 | 泄题形式 | 层级 |
| :--- | :--- | :---: |
| v1 | 前提写 `growth-arrested`，字段名含 `arrest`，与选项 `quiescence` 近义 | **词形** |
| v2 | 题面写"选择因三药共有属性而成立的因素" ⇒ 给排除法立路标 | **解法** |
| v3 | 题面用四条事实陈述疾病生物学 ⇒ 给出推理所需的全部前提 | **知识** |

三次都是"我把本该由作答者完成的那一步，替他做了"。
记录在此是为了让这个模式可被识别，而不只是修掉第三个实例。

## v4 怎么改的

| 改动 | 效果 |
| :--- | :--- |
| 整块删除四条事实清单 | 屏障状态、复发起源、低氧全部交还模型 |
| 删除三个药的机制描述，只留药名 + 数据集标签 | 复制依赖性交还模型 |
| 两个抗 VEGF 字段的选项缩为**裸方向**（`increases`/`decreases`/`unchanged`） | 原选项名含 `via_barrier_restoration`、`unchanged_because_margin_barrier_already_intact` —— **后者本身就陈述了边缘屏障完整**，又是一处泄题 |
| 新增 `notes.specification_principle` | 把界线写成显式规格，防第四次 |
| `known_gaps` 新增一条 | 诚实记录 v4 引入的新失败模式（见下） |

题面长度：**2799 字符 → 2186 字符。**

## v4 新引入的代价（已记入 `known_gaps`）

不再陈述药物机制 ⇒ **不知道 methotrexate 是 S 期特异、或不知道 doxorubicin 致死需要复制的模型，
会因为缺标准药理学而失败，而不是因为缺预期的洞察。**
选这三个药的部分理由正是这套药理学并不冷僻，但这个失败模式是真实的，读实测结果时必须考虑。

## 怎么并排对照

```bash
cd /home/gao/projects_2026Q4/48_Saul_LatchBio
A=superseded_for_comparison/evals/v3_20261105_premises_stated_the_biology

# 1) 题面逐行对照（v3 vs v4）—— 最直观，能看到整块事实清单消失
diff <(python3 -I -c "import json;print(json.load(open('$A/glioma_margin_rate_limit_v3.json'))['task'])") \
     <(python3 -I -c "import json;print(json.load(open('evals/glioma_margin_rate_limit_v4.json'))['task'])")

# 2) 题面长度对比
python3 -I -c "
import json
for p in ['$A/glioma_margin_rate_limit_v3.json','evals/glioma_margin_rate_limit_v4.json']:
    d=json.load(open(p)); print('%-40s %5d 字符'%(d['id'],len(d['task'])))
"

# 3) 读旧版可读渲染（中/英各一份；新版渲染在 tmp/）
cat "$A/rendered_v3_zh.md"
```

## 复查判据

对照时该确认的是：**v4 题面里是否还残留任何一条疾病生物学。**

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

2026-11-05 实测：以上全部为 0。`enhanc` 有 2 处命中，均在 `relapse_pattern_after_antivegf`
的选项名内（`enhancing_and_local` / `non_enhancing_and_distant`），**是答案选项而非前提，必须保留。**

## 尚未解决

**v4 同样未在任何模型上实测过。** v1→v2→v3→v4 四代缺陷全部是自审或用户追问发现的，
**没有一条来自实测**。这意味着当前对 v4 难度的全部主张都仍是判断而非证据。

归档内容一律不进客户/交付目录。
