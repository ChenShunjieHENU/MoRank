source("adaptive.R")
period <- commandArgs(TRUE)[1]
worker <- as.integer(commandArgs(TRUE)[2])
workers <- as.integer(commandArgs(TRUE)[3])
choices <- fread("settings/frozen_choices.csv")
stopifnot(nrow(choices)==2L)
new_models <- lapply(seq_len(nrow(choices)),function(i)
  read_frozen_model(file.path("models",paste0(choices$model[i],".json"))))
old_models <- lapply(choices$family,frozen_old)

score_new <- function(x) {
  values <- lapply(seq_len(nrow(choices)),function(i) {
    old <- old_models[[i]]
    new <- predict_features(new_models[[i]],x,features_for(choices$family[i],colnames(x)))
    combination(new,predict_features(old$model,x,old$cols),choices$weight[i])
  })
  setNames(as.data.frame(values),choices$reported_method)
}
top_hits <- function(score,truth,k=35L) {
  threshold <- sort(score,decreasing=TRUE)[k]
  above <- score>threshold
  ties <- score==threshold
  sum(above & truth)+(k-sum(above))*sum(ties & truth)/sum(ties)
}
metrics <- function(scores, truth, suite, design, mechanism, seed) {
  stopifnot(nrow(scores)==length(truth),all(is.finite(as.matrix(scores))))
  data.table(suite=suite,design=design,mechanism=mechanism,seed=seed,
             method=colnames(scores),auc=vapply(scores,auc,0,truth=truth),
             hits35=vapply(scores,top_hits,0,truth=truth))
}

if (period %in% c("test2","neighbor_final")) {
  selected <- fread("settings/frozen_adaptive.csv")
  specialists <- rbindlist(lapply(choices$family,function(family)
    fread(file.path("settings",paste0("frozen_specialist_",family,".csv")))))
  specialist_models <- lapply(specialists$model,function(m) {
    if (is.na(m)) return(NULL)
    read_frozen_model(file.path("models",paste0(m,".json")))
  })
  gate <- read_frozen_model("models/pool_gate.json")
  ablation_dir <- "../ablations/models"
  ablations <- readRDS(file.path(ablation_dir,"ablation_selection_frozen.rds"))
  ablation_models <- lapply(ablations$model,function(m)
    read_frozen_model(file.path(ablation_dir,paste0(m,".json"))))
  ablation_cols <- lapply(ablations$model,function(m)
    readRDS(file.path(ablation_dir,paste0(m,"_spec.rds")))$features)
  score_final <- function(f,counts) {
    prob <- predict(gate,xgb.DMatrix(matrix(pool_profile(counts),nrow=1),nthread=1))
    outputs <- list()
    for(i in seq_len(nrow(choices))) {
      family <- choices$family[i]
      pair <- predict_new_old(i,f)
      setting <- selected[family==choices$family[i]]
      weight <- with(setting,mix_weight(prob,w_nb,w_pln,rule))
      adaptive <- combination(pair$new,pair$old,weight)
      final <- pair$new
      if(prob<.5) {
        branch <- specialists[family==choices$family[i]]
        final <- pair$old
        if(branch$weight>0) {
          model <- specialist_models[[i]]
          new <- predict_features(model,f,features_for(family,colnames(f)))
          final <- combination(new,pair$old,branch$weight)
        }
      }
      outputs[[setting$reported_method]] <- adaptive
      name <- if(family=="Classical") "Classical final" else "MoRank final"
      outputs[[name]] <- final
    }
    for(i in seq_len(nrow(ablations))) {
      outputs[[ablations$family[i]]] <- predict_features(ablation_models[[i]],f,ablation_cols[[i]])
    }
    as.data.frame(outputs,check.names=FALSE)
  }
}
if (period %in% c("test","test2")) {
  # Re-run the existing count-based competitors on the same fresh count matrices.
  source("../reference/count_baselines.R")
  jobs <- jobs_for(period)
  dir.create(file.path(Sys.getenv("MORANK_AUDIT_RESULTS"),"fresh_cache"),showWarnings=FALSE,recursive=TRUE)
  ids <- seq(worker,nrow(jobs),by=workers)
  out <- lapply(ids,function(j) {
    job <- jobs[j]
    path <- file.path(Sys.getenv("MORANK_AUDIT_RESULTS"),"fresh_cache",paste0(job$seed,".rds"))
    sim <- with(job,simulate_case(suite,design,mechanism,seed))
    if (file.exists(path)) {
      record <- readRDS(path)
      f <- record$features
      count <- record$count
    } else {
      f <- both_features(sim$counts,sim$n1,sim$n0)
      count <- count_scores(sim$counts,sim$n1,sim$n0)
      saveRDS(list(features=f,count=count),path,compress=FALSE)
    }
    old <- setNames(as.data.frame(lapply(old_models,function(m)predict_features(m$model,f,m$cols))),
                    c("ClassicalStack","WSKR_A"))
    extra <- if(period=="test2") score_final(f,sim$counts) else NULL
    z <- metrics(cbind(as.data.frame(f[,baseline_columns,drop=FALSE]),old,
                       score_new(f),as.data.frame(count),extra),sim$truth,
                 job$suite,job$design,job$mechanism,job$seed)
    cat("TEST",job$suite,job$seed,"\n");flush.console()
    z
  })
}
if (period %in% c("neighbor","neighbor_final")) {
  out <- lapply(seq(worker,24L,by=workers),function(id) {
    root <- Sys.getenv("MORANK_COMMON_DIR")
    meta <- fread(file.path(root,"data","neighbor",sprintf("scores_%02d.csv",id)))
    # Scores in this comparison use the same 400 selected genes and common CPM denominator.
    x <- as.matrix(fread(file.path(root,"data","neighbor",sprintf("expression_%02d.csv.gz",id))))
    f <- cbind(distribution_features(x,51:100,1:50),local_features(x,51:100,1:50))
    mmd <- fread(file.path(root,"results",sprintf("mmdagg_neighbor_%02d.csv",id)))
    waddr <- fread(file.path(root,"results",sprintf("waddr_neighbor_%02d.csv",id)))
    stopifnot(identical(meta$gene_index,mmd$gene_index),identical(meta$gene_index,waddr$gene_index))
    scores <- cbind(as.data.frame(meta[,10:ncol(meta),with=FALSE]),score_new(f),
                    MMDAgg=-mmd$inverted_level,waddR=ifelse(is.na(waddr$pval),-1,-waddr$pval))
    if (period=="neighbor_final") {
      sim <- simulate_case("NB",1L,match(meta$scenario[1],mechanisms),meta$seed[1])
      full_x <- log2(sweep(sim$counts,2,colSums(sim$counts),"/")*1e6+1)
      stopifnot(max(abs(full_x[meta$gene_index,]-x))<1e-11,
                identical(sim$truth[meta$gene_index],meta$truth))
      full_f <- cbind(distribution_features(full_x,51:100,1:50),local_features(full_x,51:100,1:50))
      final <- score_final(full_f,sim$counts)
      scores <- cbind(scores,final[meta$gene_index,,drop=FALSE])
    }
    metrics(scores,meta$truth,"Neighbor",1L,match(meta$scenario[1],mechanisms),meta$seed[1])
  })
}
fwrite(rbindlist(out),file.path(Sys.getenv("MORANK_AUDIT_RESULTS"),paste0(period,"_",worker,".csv")))
cat("DONE",period,worker,length(out),"datasets\n")
