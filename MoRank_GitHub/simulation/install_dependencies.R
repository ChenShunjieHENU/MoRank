cran <- "https://cloud.r-project.org"
packages <- c("data.table", "xgboost", "matrixStats", "Matrix", "jsonlite",
              "igraph", "RSpectra", "cli", "R.utils", "statmod", "BiocManager")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = cran)
bioc <- c("DESeq2", "edgeR", "waddR")
missing <- setdiff(bioc, rownames(installed.packages()))
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
