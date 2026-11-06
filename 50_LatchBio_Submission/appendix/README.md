# Appendix — LLM testing

## 要做什么（约 10 分钟）

把 `prompt_problem1.txt` 与 `prompt_problem2.txt` 的**全文**分别粘进 LLM，每题至少测 2 次。
免费账号即可。建议用两个不同产品（例如 ChatGPT 免费版 + Gemini 免费版），
这样"哪里做得好、哪里推理失败"的观察更有说服力。

**每次新开一个空白对话**，不要在同一会话里连着问两题（前一题会污染后一题）。

## 要记录什么

每次测试存成一个文件：

| 文件名 | 内容 |
| :--- | :--- |
| `p1_<产品名>_run1.txt` | 第一行写 `MODEL: <产品与版本>  DATE: 2026-11-05`，其后粘贴**完整响应**，一个字不删 |

完整响应是硬性要求，不是摘要。

## 怎么判分

```bash
cd /home/gao/projects_2026Q4/50_LatchBio_Submission
python3 -I scripts/grade.py p1 < appendix/p1_chatgpt_run1.txt
python3 -I scripts/grade.py p2 < appendix/p2_gemini_run1.txt
```

输出 `PASS` / `FAIL` 加失败原因。P2 若模型答 `unchanged`，脚本会打印 `[REVIEW]`
提示需人工看 rationale —— 这是刻意的，那一条无法机械判定。

## 文件清单

- `prompt_problem1.txt` — 给模型的确切 prompt（与正文 §3 逐字相同）
- `prompt_problem2.txt` — 给模型的确切 prompt（与正文 §4 逐字相同）
- `p1_*.txt` / `p2_*.txt` — 完整响应，待你粘贴
