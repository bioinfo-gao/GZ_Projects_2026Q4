# _common.R —— 每一章开头 source() 一次
# 统一：数据目录、下载函数、随机种子、ggplot2 主题

options(
  timeout = max(1800, getOption("timeout")),   # 大文件下载
  width = 90,
  digits = 4,
  dplyr.summarise.inform = FALSE
)

DATA_DIR <- "data"
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

#' 只在本地不存在时下载；返回本地路径
fetch <- function(url, dest = file.path(DATA_DIR, basename(url))) {
  if (!file.exists(dest)) {
    message("Downloading ", url)
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

#' 缓存昂贵计算：表达式只在文件不存在时求值
cached <- function(path, expr) {
  path <- file.path(DATA_DIR, path)
  if (file.exists(path)) return(readRDS(path))
  value <- expr            # 懒求值：只有走到这里才会真的计算
  saveRDS(value, path)
  value
}

set.seed(2026)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 12))
}
