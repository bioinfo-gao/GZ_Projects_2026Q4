# v1_20261126_knn_nondeterministic（P47 round-1 描述性，第三次链运行产物，11-26 00:04–00:59）

**缺陷**
1. step14 kNN 图用 pynndescent（scanpy 默认近似近邻，多线程）——同 random_state 重建图不可复现
   （实测 nnz 2590492 vs 2590612；leiden 簇数 10 / 9 / 11，ARI 0.71–0.86）。leiden 在固定图上本身确定（重载后两遍 ARI=1.0）。
   ⇒ 本版 leiden 簇号、簇组成表、15 decontX(z=leiden)、16 sample×cluster pseudobulk 均依赖一次不可复现的图。
2. step18 伪重复自证断言前提写错：同一样本随机对半得 0 个显著基因是**正确的校准结果**，不是伪重复的证据。
3. step19 `assert_no_palette_collision` 传入 dict（迭代出键名）→ 崩溃，未出任何图。

**修法**：step14 `sc.pp.neighbors(..., transformer="sklearn")` 精确 kNN；step18 确定性检查改为从 ckpt 嵌入重建图并逐元素比对；伪重复改为同名义条件两个样本的细胞级检验，随机对半保留为校准；step19 传 `.values()`。

**对照**：`diff scripts/<f> ../../../scripts/<f>`；产物与 `output_results/round1_descriptive_1125/` 同名文件并排比。
旧 ckpt 原样改名保留：`ckpt_1125/ckpt_14_clustered_v1_20261126_knn_nondeterministic.h5ad`（原始文件本身，非复现）。
样本级 pseudobulk（16_pseudobulk_sample_counts）与 17 对比不依赖 leiden，预期新旧一致——可作为回归检查。
