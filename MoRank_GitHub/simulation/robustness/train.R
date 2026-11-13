suppressPackageStartupMessages({library(data.table);library(xgboost);library(jsonlite)})
setDTthreads(1)
source("features.R")
args <- commandArgs(TRUE)
mode <- args[1]
features_for <- function(family, all_names) {
  switch(family, Classical=classic_features, MoRank37=all_names[1:37], MoRank48=all_names)
}
load_split <- function(split) {
  jobs <- jobs_for(split)
  units <- lapply(jobs$seed,function(seed) readRDS(file.path(Sys.getenv("MORANK_MIXED_CACHE","cache"),split,paste0(seed,".rds"))))
  stopifnot(identical(vapply(units,function(x) x$job$seed,0),as.numeric(jobs$seed)))
  units
}
frozen_old <- function(family) {
  name <- if(family=="Classical") "ClassicalStack_d4_t400" else "WSKR_A_d4_t400"
  dir <- "../reference/models"
  list(model=read_frozen_model(file.path(dir,paste0(name,".json"))),
       cols=readRDS(file.path(dir,paste0(name,"_spec.rds")))$features)
}
predict_features <- function(model,x,cols) {
  predict(model,xgb.DMatrix(x[,cols,drop=FALSE],nthread=1))
}
rank_fraction <- function(score) (rank(score,ties.method="average")-.5)/length(score)
combination <- function(new, old, weight) {
  if(weight==1) return(new)
  weight*rank_fraction(new)+(1-weight)*rank_fraction(old)
}
scores_by_unit <- function(scores, units) {
  aucs <- vapply(seq_along(units),function(i) auc(scores[[i]],units[[i]]$truth),0)
  info <- rbindlist(lapply(units,function(x)x$job))
  info[,auc:=aucs]
  info
}

