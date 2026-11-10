# P52 research plan and review

Plan Date（创建）: 2026-11-09

## 更新记录
- 2026-11-09 — 首次清点；尚未正式运行，等待科学口径确认。

## 目标和输入
用户授权全过程：FASTQ → QC → nf-core/rnaseq → DESeq2 → GO/KEGG ORA + GSEA → 客户报告与交付包（仅准备，待用户审阅）。
项目: /home/gao/projects_2026Q4/52_Lijian_Wu_RNAseq_by_Codex
原始目录: /home/gao/Dropbox/Quote_09292601
截至 2026-11-09 本轮清点，该目录可见 48 个 FASTQ、24 对，压缩总量 135.823 GiB。文件完整性及 read pairs 尚未验证；文件存在不等于该批同步完毕。
样本与每样本数据量 SSOT: ../scripts/0_sample_table.csv。
送样单确认 Mouse、8 组各 3 样本，并列 6 个比较。

## 正式运行前待确认
1. FASTQ Control/Ad/PGE/TNFA 对应表单 un/adenosine/PEG2（比较栏 pge2）/TNF，重复编号一一对应。
2. 表单 log2FC cutoff 写 Use default (1.0)，与本机 skill 默认 0.263 不同；需选择本项目正式口径。
3. 多组过滤建议：本对比两组 ≥3 个样本 estimated counts ≥10，或全体 ≥30% 样本（24 样本取整为 ≥8）≥10；每对比独立检验集合，independentFiltering=FALSE，集合内 BH。

## Evidence decision table
| 研究问题 | 观察/证据单位 | 保留与去重规则及理由 | 判定指标 | 允许的结论 |
| :--- | :---: | :---: | :---: | :---: |
| 组间表达差异 | 生物学样本；基因层 Salmon estimated counts | 全部24样本；不按read数当生物学重复；不因PCA位置自动剔除 | 六项预设比较，DESeq2，BH ≤0.05，fold-change待确认 | 组间表达关联，不推断因果 |
| 低表达控制 | 每对比基因集合 | 上述多组过滤待确认；删除该规则会增加低信息检验，过严过滤会漏掉处理特异表达 | 保存每基因通过/排除原因及前后数量 | 未通过过滤不等于无表达 |
| 比对及定量 | paired fragment与alignment不同层次 | 保留pipeline alignment与QC记录；Salmon处理多位置归属，不将secondary当独立分子 | mapping、strand、assigned counts及样本质量审计 | 定量质量与定位局限 |
| 通路关联 | 背景为该对比可检验且可映射基因 | ORA分方向、GSEA用完整有序可检验集合；审计ID丢失 | BH校正与富集方向 | 通路关联，不等于活性实测 |

## 执行与完成标准
确认后编写编号脚本并检查；使用独立新目录和tmux、资源总上限28核/56线程，先核查其他任务。正式启动前配置watchdog与中心active登记。当前未启动正式程序，不宣称分析完成。完成要求含表达矩阵、QC、六比较、富集、PDF图与点数核验、英文报告、主张到数据来源对照及交付清单。
