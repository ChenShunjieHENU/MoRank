suppressPackageStartupMessages({library(data.table);library(xgboost);library(jsonlite)})
setDTthreads(1)
source("features.R")
args <- commandArgs(TRUE)
mode <- args[1]
audit_path <- Sys.getenv("MORANK_FIT_OUTPUT")
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

if (mode=="train") {
  family <- args[2]
  fit <- load_split("fit")
  tune <- load_split("tune")
  cols <- features_for(family,colnames(fit[[1]]$features))
  x <- do.call(rbind,lapply(fit,function(x)x$features[,cols,drop=FALSE]))
  y <- as.integer(unlist(lapply(fit,function(x)x$truth)))
  dm <- xgb.DMatrix(x,label=y,nthread=6)
  rm(fit,x,y);gc()
  previous <- frozen_old(family)
  old_scores <- lapply(tune,function(u)predict_features(previous$model,u$features,previous$cols))
  tune_mats <- lapply(tune,function(u)xgb.DMatrix(u$features[,cols,drop=FALSE],nthread=1))
  grid <- CJ(depth=c(3L,5L),rounds=c(160L,360L))
  options <- list()
  for (i in seq_len(nrow(grid))) {
    depth <- grid$depth[i]; rounds <- grid$rounds[i]
    params <- list(objective="binary:logistic",tree_method="hist",eval_metric="logloss",
       nthread=6,max_depth=depth,eta=.05,min_child_weight=20,lambda=10,
       subsample=.8,colsample_bytree=.9,seed=812345L)
    model <- xgb.train(params,dm,nrounds=rounds,verbose=0)
    name <- sprintf("%s_d%d_t%d",family,depth,rounds)
    scores <- lapply(tune_mats,function(d)predict(model,d))
    xgb.save(model,file.path(audit_path,"models",paste0(name,".json")))
    for (weight in c(.5,.75,1)) {
      score <- Map(function(new,old)combination(new,old,weight),scores,old_scores)
      cells <- scores_by_unit(score,tune)[,.(auc=mean(auc)),by=.(suite,design,mechanism)]
      z <- dcast(cells[,.(auc=mean(auc)),by=suite],. ~ suite,value.var="auc")
      z[["."]] <- NULL
      z[, ':='(family=family,model=name,depth=depth,rounds=rounds,weight=weight)]
      z[,selection_auc:=mean(unlist(.SD)),.SDcols=suites]
      options[[length(options)+1L]] <- z
    }
    fwrite(rbindlist(options),file.path(audit_path,"results",paste0("tune_",family,".csv")))
    cat("TUNED",name,round(max(rbindlist(tail(options,3))$selection_auc),5),"\n")
    flush.console()
  }
}
if (mode=="select") {
  files <- file.path(audit_path,"results",paste0("tune_",c("Classical","MoRank37","MoRank48"),".csv"))
  grid <- rbindlist(lapply(files,fread))
  stopifnot(nrow(grid)==36L)
  pick <- function(z) {
    close <- z[selection_auc >= max(selection_auc)-.0005]
    close[,features:=ifelse(family=="MoRank48",48L,ifelse(family=="MoRank37",37L,9L))]
    setorder(close,features,depth,rounds,-weight,-selection_auc)
    close[1]
  }
  chosen <- rbind(pick(grid[family=="Classical"]),pick(grid[family!="Classical"]))
  chosen[,reported_method:=c("Classical retrained","MoRank optimized")]
  fwrite(chosen,file.path(audit_path,"results/frozen_choices.csv"))
  print(chosen)
}

