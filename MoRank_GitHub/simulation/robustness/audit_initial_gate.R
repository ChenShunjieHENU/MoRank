source("train.R")
mode <- commandArgs(TRUE)[1]
dest <- Sys.getenv("MORANK_AUDIT_GATE_OUTPUT")
choices <- fread("settings/frozen_choices.csv")
old_models <- lapply(choices$family,frozen_old)
new_models <- lapply(choices$model,function(model)
  read_frozen_model(file.path("models",paste0(model,".json"))))

pool_profile <- function(counts) {
  cpm <- log2(sweep(counts,2,pmax(colSums(counts),1),"/")*1e6+1)
  mu <- rowMeans(cpm)
  sd <- matrixStats::rowSds(cpm)
  zero <- rowMeans(counts==0)
  cv <- sd / pmax(mu,.01)
  distribution <- function(x) as.numeric(quantile(x,c(.1,.25,.5,.75,.9),names=FALSE))
  c(distribution(mu),distribution(sd),distribution(zero),distribution(cv),
    mean(zero),mean(sd),median(sd[mu>=median(mu)]))
}
profile_matrix <- function(split) {
  jobs <- jobs_for(split)
  do.call(rbind,lapply(seq_len(nrow(jobs)),function(i) {
    sim <- with(jobs[i],simulate_case(suite,design,mechanism,seed,split!="test"))
    pool_profile(sim$counts)
  }))
}
predict_new_old <- function(i,x) {
  old <- old_models[[i]]
  list(new=predict_features(new_models[[i]],x,features_for(choices$family[i],colnames(x))),
       old=predict_features(old$model,x,old$cols))
}
mix_weight <- function(p,w_nb,w_pln,rule) {
  if(rule=="Hard") p <- as.numeric(p>=.5)
  w_nb + p*(w_pln-w_nb)
}

if(mode=="train_gate") {
  dfit <- profile_matrix("fit")
  dval <- profile_matrix("tune")
  jobs <- jobs_for("fit")
  targets <- startsWith(jobs$suite,"PLN")
  gate <- xgb.train(list(objective="binary:logistic",tree_method="hist",max_depth=2,eta=.07,
                          min_child_weight=5,lambda=10,nthread=1,seed=832345L),
                     xgb.DMatrix(dfit,label=targets),nrounds=120L)
  xgb.save(gate,file.path(dest,"models/pool_gate.json"))
  prob <- predict(gate,xgb.DMatrix(dval))
  tune <- load_split("tune")
  info <- jobs_for("tune")[,.(suite,design,mechanism,seed)]
  info[,pln_probability:=prob]
  fwrite(info,file.path(dest,"results/tune_domain_weights.csv"))
  macros <- function(score) {
    z <- scores_by_unit(score,tune)[,.(auc=mean(auc)),by=.(suite,design,mechanism)]
    v <- z[,.(auc=mean(auc)),by=suite]
    setNames(v$auc,v$suite)[suites]
  }
  options <- list()
  for(i in seq_len(nrow(choices))) {
    pairs <- lapply(tune,function(u)predict_new_old(i,u$features))
    new <- lapply(pairs,function(x)x$new)
    old <- lapply(pairs,function(x)x$old)
    target <- pmax(macros(new),macros(old))
    grid <- CJ(w_nb=c(0,.25,.5,.75,1),w_pln=c(.75,1),rule=c("Soft","Hard"))
    for(j in seq_len(nrow(grid))) {
      w <- with(grid[j],mix_weight(prob,w_nb,w_pln,rule))
      result <- macros(Map(combination,new,old,w))
      options[[length(options)+1L]] <- cbind(
        data.table(family=choices$family[i]),grid[j],
        as.data.table(as.list(result)),
        data.table(min_gain=min(result-target),avg_auc=mean(result)))
    }
  }
  all <- rbindlist(options)
  fwrite(all,file.path(dest,"results/adaptive_tune.csv"))
  setorder(all,family,-min_gain,-avg_auc,rule)
  selected <- all[, .SD[1], by=family]
  selected[,reported_method:=ifelse(family=="Classical","Classical adaptive","MoRank adaptive")]
  fwrite(selected,file.path(dest,"results/frozen_adaptive.csv"))
  print(selected)
}
