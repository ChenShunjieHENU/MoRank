source("adaptive.R")
family <- commandArgs(TRUE)[2]
audit_path <- Sys.getenv("MORANK_FIT_OUTPUT")
stopifnot(family %in% c("Classical","MoRank37"))
fit <- load_split("fit")
fit <- Filter(function(u)!startsWith(u$job$suite,"PLN"),fit)
archived <- list.files(file.path(Sys.getenv("MORANK_DEVELOPMENT"), "fit"),pattern="\\.rds$",full.names=TRUE)
stopifnot(length(archived)==144L)
fit <- c(fit,lapply(archived,readRDS))
tune <- Filter(function(u)!startsWith(u$job$suite,"PLN"),load_split("tune"))
cols <- features_for(family,colnames(fit[[1]]$features))
x <- do.call(rbind,lapply(fit,function(u)u$features[,cols,drop=FALSE]))
y <- as.integer(unlist(lapply(fit,function(u)u$truth)))
dm <- xgb.DMatrix(x,label=y,nthread=8)
rm(x,y,fit);gc()
old <- frozen_old(family)
old_scores <- lapply(tune,function(u)predict_features(old$model,u$features,old$cols))
mats <- lapply(tune,function(u)xgb.DMatrix(u$features[,cols,drop=FALSE],nthread=1))
macro <- function(score) {
  a <- scores_by_unit(score,tune)[,.(auc=mean(auc)),by=.(suite,design,mechanism)]
  a <- a[,.(auc=mean(auc)),by=suite]
  setNames(a$auc,a$suite)[c("NB","Sensitivity")]
}
reference <- macro(old_scores)
options <- list(data.table(model=NA_character_,depth=0L,rounds=0L,weight=0,
                           NB=reference[1],Sensitivity=reference[2],min_gain=0,avg_gain=0))
grid <- CJ(depth=c(4L,6L),rounds=c(400L,750L))
for(i in seq_len(nrow(grid))) {
  d <- grid$depth[i]; t <- grid$rounds[i]
  model <- xgb.train(list(objective="binary:logistic",tree_method="hist",eval_metric="logloss",
    nthread=8,max_depth=d,eta=.05,min_child_weight=20,lambda=10,subsample=.8,
    colsample_bytree=.9,seed=854321L),dm,nrounds=t,verbose=0)
  name <- sprintf("%s_NB_d%d_t%d",family,d,t)
  xgb.save(model,file.path(audit_path,"models",paste0(name,".json")))
  scores <- lapply(mats,function(a)predict(model,a))
  for(w in c(.25,.5,.75,1)) {
    result <- macro(Map(function(new,previous)combination(new,previous,w),scores,old_scores))
    options[[length(options)+1L]] <- data.table(
      model=name,depth=d,rounds=t,weight=w,NB=result[1],Sensitivity=result[2],
      min_gain=min(result-reference),avg_gain=mean(result-reference))
  }
  cat("TUNED",name,"\n");flush.console()
}
all <- rbindlist(options)
all[,family:=family]
fwrite(all,file.path(audit_path,"results",paste0("specialist_grid_",family,".csv")))
setorder(all,-min_gain,-avg_gain,depth,rounds)
fwrite(all[1],file.path(audit_path,"results",paste0("frozen_specialist_",family,".csv")))
print(all[1:4])

