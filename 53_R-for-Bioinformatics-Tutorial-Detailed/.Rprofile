# p53 项目级 .Rprofile —— 补上 conda 激活时才有的 4 个环境变量
# 为什么：Positron 直接启动 p53_rbio/bin/R、不做 conda activate ⇒ `conda env config vars` 不生效。
# 只在变量缺失时才设（终端里 mamba activate 过就原样保留）。数值须与 install-conda.sh 中
# `conda env config vars set` 那一行一致（那里是源头；改一处必须同改另一处）。
# 仅当 R 启动时工作目录 = 本项目根才会被读取；说明见 ~/Utilities/IDE/Positron/R_setting_in_server/Rprofile_for_conda_R_in_Positron.md
local({
  want <- c(
    OPENBLAS_NUM_THREADS = "8",
    OMP_NUM_THREADS      = "8",
    R_USER_CACHE_DIR     = "/mnt/ex_8T_SSD/R_user_cache",
    XDG_CACHE_HOME       = "/Work_bio/gao/configs/.cache"
  )
  missing <- want[Sys.getenv(names(want)) == ""]
  if (length(missing)) do.call(Sys.setenv, as.list(missing))

  # 线程变量必须在 BLAS 载入前设置才有效；R 启动时 OpenBLAS 已载入 ⇒ 运行中再限一次
  n <- as.integer(Sys.getenv("OPENBLAS_NUM_THREADS"))
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(n)
    RhpcBLASctl::omp_set_num_threads(as.integer(Sys.getenv("OMP_NUM_THREADS")))
  }
  if (interactive())
    message(sprintf("[p53 .Rprofile] 补设 %d 个变量 (%s); BLAS 线程=%s",
                    length(missing), paste(names(missing), collapse = ","),
                    if (requireNamespace("RhpcBLASctl", quietly = TRUE))
                      RhpcBLASctl::blas_get_num_procs() else "未限制(缺 RhpcBLASctl)"))
})

# 项目级 .Rprofile 存在时 R 不再读 ~/.Rprofile；若以后建了用户级的，在这里显式读入
if (file.exists("~/.Rprofile")) source("~/.Rprofile")
