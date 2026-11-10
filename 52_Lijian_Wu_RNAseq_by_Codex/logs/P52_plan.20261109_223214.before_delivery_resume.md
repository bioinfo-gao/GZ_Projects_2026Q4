# P52 research plan and review

Plan Date（创建）: 2026-11-09

## 更新记录
- 2026-11-09 — 首次清点；当时等待科学口径确认。
- 2026-11-09 — 用户确认样本映射、|log2FC| ≥0.263及大样本过滤规则；RNA-seq skill中心源文件已同步更新。

## 目标和输入
用户授权全过程：FASTQ → QC → nf-core/rnaseq → DESeq2 → GO/KEGG ORA + GSEA → 客户报告与交付包（仅准备，待用户审阅）。
项目: /home/gao/projects_2026Q4/52_Lijian_Wu_RNAseq_by_Codex
原始目录: /home/gao/Dropbox/Quote_09292601
截至 2026-11-09 本轮清点，该目录可见 48 个 FASTQ、24 对，压缩总量 135.823 GiB。文件完整性及 read pairs 尚未验证；文件存在不等于该批同步完毕。
样本与每样本数据量 SSOT: ../scripts/0_sample_table.csv。
送样单确认 Mouse、8 组各 3 样本，并列 6 个比较。

## 已确认设置（用户确认：2026-11-09）
1. FASTQ Control/Ad/PGE/TNFA 对应表单 un/adenosine/PEG2（比较栏 pge2）/TNF，重复编号一一对应；用户已确认。
2. 表单 log2FC cutoff 写 Use default (1.0)，与本机 skill 默认 0.263 不同；本项目由用户明确指定采用 |log2FC| ≥0.263（覆盖表单的1.0）。
3. 多组过滤建议：本对比两组 ≥3 个样本 estimated counts ≥10，或全体 ≥1/3 样本（24 样本为 ≥8）≥10；每对比独立检验集合，independentFiltering=FALSE，集合内 BH。

### 通用人数公式
总样本数 N>12：`sum(counts>=10 across all samples)>=ceil(N/3) OR sum(counts>=10 across contrast samples)>=ceil(n_contrast/2)`。表达阈值为≥10 estimated counts；人数向上取整。本项目N=24、n_contrast=6，即≥8或≥3。

### 六项预设比较（左侧为实验组，右侧为基线）
- Crtc2 vs Control
- Crtc2_TNFA vs Control_TNFA
- Crtc2_PGE vs Control_PGE
- Crtc2_Ad vs Control_Ad
- Control_PGE vs Control
- Control_Ad vs Control

## Evidence decision table
| 研究问题 | 观察/证据单位 | 保留与去重规则及理由 | 判定指标 | 允许的结论 |
| :--- | :---: | :---: | :---: | :---: |
| 组间表达差异 | 生物学样本；基因层 Salmon estimated counts | 全部24样本；不按read数当生物学重复；不因PCA位置自动剔除 | 六项预设比较，DESeq2，BH ≤0.05，|log2FC| ≥0.263 | 组间表达关联，不推断因果 |
| 低表达控制 | 每对比基因集合 | 上述多组过滤已由用户确认；删除该规则会增加低信息检验，过严过滤会漏掉处理特异表达 | 保存每基因通过/排除原因及前后数量 | 未通过过滤不等于无表达 |
| 比对及定量 | paired fragment与alignment不同层次 | 保留pipeline alignment与QC记录；Salmon处理多位置归属，不将secondary当独立分子 | mapping、strand、assigned counts及样本质量审计 | 定量质量与定位局限 |
| 通路关联 | 背景为该对比可检验且可映射基因 | ORA分方向、GSEA用完整有序可检验集合；审计ID丢失 | BH校正与富集方向 | 通路关联，不等于活性实测 |

## 执行与完成标准
已获用户全过程授权及科学口径确认，接下来编写编号脚本并检查；使用独立新目录和tmux、资源总上限28核/56线程，先核查其他任务。正式启动前配置watchdog与中心active登记。当前未启动正式程序，不宣称分析完成。完成要求含表达矩阵、QC、六比较、富集、PDF图与点数核验、英文报告、主张到数据来源对照及交付清单。

## 执行状态（2026-11-09 10:35 EST）
- 已启动 tmux `p52_full`：8路 gzip -t → nf-core/rnaseq 3.15.1 → FastQC R1/R2计数一致性 → DESeq2六比较 → GO/KEGG ORA + GSEA；R语法与依赖加载检查已通过。
- 总控日志 `../logs/master.log`，状态 `../logs/MASTER_STATUS`；监控 `p52_watch`，每10分钟，cron每30分钟复活；中心active与Codex本会话续跑已登记。
- 资源：CPU线程总预算56、内存108GB，STAR每任务18线程/32GB、最多3并发；启动实测available约113GiB，项目盘剩余约2TiB。预计数小时至一天，以首批实测回填。
- 新目录 output_results、work、Data_Analysis_YYYYMMDD；不改原始FASTQ、不复用P51结果。复用P51脚本结构，已修改为本项目输入、0.263及1/3 OR 1/2过滤；所有基因及每比较过滤原因写machine-readable audit。
- 监控曾将尚无Nextflow日志的FASTQ预检阶段误报STALLED，已记录并修正项目薄封装：跨阶段禁用不适用的Nextflow心跳检查，保留session/资源/long-task检查并显式指定canonical mailer。Nextflow日志写项目根。LONG_TASK模拟7小时任务实测exit15；此前日志保留。
- 正式结果尚未生成；不会使用P51报告的生物学结论。结果生成后解读、图形审阅、报告、交付审计仍为本全过程的必需步骤。
