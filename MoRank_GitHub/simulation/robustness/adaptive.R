source("train.R")
mode <- commandArgs(TRUE)[1]
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

