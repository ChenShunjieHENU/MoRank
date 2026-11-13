suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
source("features.R")
args <- commandArgs(TRUE)
split <- args[1]
if (split=="test2") source("../reference/count_baselines.R")
worker <- as.integer(args[2])
workers <- as.integer(args[3])
jobs <- jobs_for(split)
stopifnot(!anyDuplicated(jobs$seed))
dir.create(file.path(Sys.getenv("MORANK_MIXED_CACHE","cache"),split),recursive=TRUE,showWarnings=FALSE)
for (j in seq(worker,nrow(jobs),by=workers)) {
  job <- jobs[j]
  path <- file.path(Sys.getenv("MORANK_MIXED_CACHE","cache"),split,paste0(job$seed,".rds"))
  if (file.exists(path)) next
  sim <- with(job,simulate_case(suite,design,mechanism,seed,split %in% c("fit","tune")))
  f <- both_features(sim$counts,sim$n1,sim$n0)
  stopifnot(nrow(f)==length(sim$truth),ncol(f)==48L)
  # Keep the validation scores and labels in distinct dataset records.
  record <- list(job=job,features=f,truth=sim$truth)
  if (split=="test2") record$count <- count_scores(sim$counts,sim$n1,sim$n0)
  saveRDS(record,path,compress=FALSE)
}
cat("DONE",split,worker,"\n")
