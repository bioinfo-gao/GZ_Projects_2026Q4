# P46 — 致 Yeji Park：开跑前两个必答问题（2026-11-06 草稿）

- **状态**：**已 APPEND 进 Gmail Drafts，未发送**（实测 `before=0 after=1 uid=1`）。
  等用户逐字过目后由用户本人发出。
- **起草脚本（SSOT，正文只存在于这一处）**：`scripts/1_append_client_email_1106.py`
  —— 正文与硬断言都在脚本里，本文不复制正文，避免两处漂移。
- **机械化检查**：`python3 ~/.claude/hooks/checks/check_client_email.py scripts/1_append_client_email_1106.py`
  ⇒ exit 0；脚本自带 13 条断言（含 4 条正向断言）全过。
- **To**：Yeji.Park@childrens.harvard.edu
  **Cc**：j.shi@athenomics.com、kchung@athenomics.com、jianjun@athenomics.com
- **线程**：`19fdd85e1110e93e`（"Re: Inquiry About Whole Genome Sequencing Services [EXTERNAL]"，
  实测 30 封）。`In-Reply-To` = 客户 2026-11-06 18:11 UTC 那封（带样本信息文件的）
  `<SA1PR16MB6978A079A13EDD94617380AFAD762@SA1PR16MB6978.namprd16.prod.outlook.com>`，
  Message-ID 由 IMAP 实测取得，非推测。

## 信里问了什么

| # | 问题 | 对应计划 | 为什么必须问 |
| :--- | :--- | :---: | :--- |
| 1 | 本轮 Cas9 以何形式进入细胞（蛋白 + gRNA / 瞬时质粒 / 克隆自带整合的 Cas9 构建体）；若是新构建体请给图谱 | §9 第 1 项 → Stage 3 | 整合筛查**只能找到我们给了参考序列的东西**；不给参考即系统性盲区。信里已说明我方手上已有 `mCherry-BP1-2 pLPC-Puro` 与 `pMSCVpuro-Cas9-EGFP` 两份图谱、两者都会查 ⇒ 客户只需回答"是否另有构建体" |
| 2 | 样本分组确认（Ctr_1–3 / Sample 1–4 单等位 / Sample 5–7 双等位） | §9 第 2 项 → §7.2 | 问卷两套编号（1–7 与 1–10），我方读作偏移 3。信里**明确写成"we read them as"**，是我方读法待客户确认，不写成既定事实 |
| 3 | （不是问题，是期望管理）十个克隆各自单细胞扩增、各带私有变异；三个对照用于估计该本底；结论将落在"编辑克隆的新增变化是否多于对照"这一层，**不会把单个克隆的单个变化归因于编辑** | §8.4 的 decision table | 客户问卷写的是 "analyze any further mutations by the CRISPR system"。把可得结论的层级**在出结果之前**讲清楚，避免结果阶段的归因落差 |

## 刻意没写进去的内容

| 没写 | 理由 |
| :--- | :--- |
| 浓度 `total 200ng/ul` 的歧义 | 用户 2026-11-06 裁定关闭：Hung 已完成建库，问了不改变任何决定（§9 第 3 项） |
| 测序数据到没到、什么时候到 | 「不提数据传输」—— 提了等于催，而客户无法作答 |
| gRNA2 的 19 nt vs 20 nt | 两轮未答，两版都做的成本极低（P31 已有脚本），报告里并列即可 |
| 与 P31 两样本合并分析 | 用户 2026-11-06 裁定 **不做**（§9 第 4 项） |
| 任何价格、工时、免费范围 | 商务口径归 Jing；10-29 实发邮件已划界，本信不动它 |
