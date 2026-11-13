source("adaptive.R")

args <- commandArgs(TRUE)
phase <- args[1]
destination <- Sys.getenv("MORANK_AUDIT_REPLICATION")
sample_designs <- rbind(primary_designs, small_designs)
families <- c("NB", "Sensitivity", "PLN")

profile_jobs <- function(split) {
  base <- if (split == "profile_fit") 871000000L else 872000000L
  jobs <- rbindlist(lapply(families, function(family) {
    z <- as.data.table(expand.grid(design=1:6,
                  mechanism=seq_len(if (family=="PLN") 4L else 8L), replicate=1:2))
    z[, family := family]
    z
  }))
  jobs[, seed := base + match(family,families)*100000L + design*1000L + mechanism*10L + replicate]
  jobs[]
}

simulate_profile <- function(family, design, mechanism, seed) {
  n <- sample_designs[design,]
  if (family=="PLN") return(generate_pln(mechanism,n[1],n[2],seed))
  set.seed(seed)
  config <- if (family=="Sensitivity") stress_config(TRUE) else default_config()
  max_n <- max(n)
  sim <- simulate_pool_fast(mechanisms[mechanism],max_n,config)
  cols <- c(seq_len(n[2]), max_n+seq_len(n[1]))
  list(counts=sim$counts[,cols,drop=FALSE], truth=sim$truth$is_signal, n1=n[1], n0=n[2])
}

stratified_profile <- function(counts) {
  depth <- pmax(colSums(counts),1)
  depth <- depth/mean(depth)
  normalized <- sweep(counts,2,depth,"/")
  avg <- rowMeans(normalized)
  variance <- matrixStats::rowVars(normalized)
  shot <- rowMeans(sweep(counts,2,depth^2,"/"))
  excess <- log1p(pmax(variance-shot,0)/pmax(avg^2,.01))
  zero <- rowMeans(counts==0)
  x <- log1p(normalized)
  centred <- x-rowMeans(x)
  sd <- sqrt(rowMeans(centred^2))
  skew <- rowMeans((centred/pmax(sd,1e-8))^3)
  # Within-abundance summaries retain the relationship between mean and dispersion.
  bin <- pmin(5L, floor(5*(rank(avg,ties.method="average")-1)/length(avg))+1L)
  summarize <- function(value) as.numeric(quantile(value,c(.1,.5,.9),names=FALSE))
  c(pool_profile(counts), unlist(lapply(1:5,function(i) {
    use <- bin==i
    if (!any(use)) return(rep(0,9))
    c(summarize(excess[use]),summarize(zero[use]),summarize(skew[use]))
  })))
}

if (phase %in% c("profile_fit","profile_tune")) {
  worker <- as.integer(args[2])
  workers <- as.integer(args[3])
  jobs <- profile_jobs(phase)
  table <- rbindlist(lapply(seq(worker,nrow(jobs),by=workers),function(i) {
    job <- jobs[i]
    sample <- with(job,simulate_profile(family,design,mechanism,seed))
    x <- stratified_profile(sample$counts)
    stopifnot(length(x)==68L,all(is.finite(x)))
    cbind(job,as.data.table(as.list(setNames(x,paste0("x",seq_along(x))))))
  }))
  fwrite(table,file.path(destination,"results",sprintf("%s_%02d.csv",phase,worker)))
  cat(phase,worker,nrow(table),"completed\n")
}

if (phase=="select_gate") {
  read_profiles <- function(split) {
    paths <- list.files(file.path(destination,"results"),
                        pattern=paste0("^",split,"_[0-9]+\\.csv$"),full.names=TRUE)
    z <- rbindlist(lapply(paths,fread))
    stopifnot(nrow(z)==nrow(profile_jobs(split)),!anyDuplicated(z$seed),
              setequal(z$seed,profile_jobs(split)$seed))
    z
  }
  train <- read_profiles("profile_fit")
  tune <- read_profiles("profile_tune")
  stopifnot(length(intersect(train$seed,tune$seed))==0L)
  train[, is_pln := as.integer(family=="PLN")]
  tune[, is_pln := as.integer(family=="PLN")]
  train[, row_weight := .5/.N, by=is_pln]
  train[, row_weight := row_weight*nrow(train)]
  cols <- paste0("x",1:68)
  old_gate <- read_frozen_model("models/pool_gate.json")
  probability <- list(Original=predict(old_gate,xgb.DMatrix(as.matrix(tune[,cols[1:23],with=FALSE]),nthread=1)))
  for (mode in c("Balanced","Stratified")) {
    use <- if (mode=="Balanced") cols[1:23] else cols
    model <- xgb.train(list(objective="binary:logistic",tree_method="hist",max_depth=2L,
      eta=.07,min_child_weight=5,lambda=10,nthread=2,seed=871456L),
      xgb.DMatrix(as.matrix(train[,use,with=FALSE]),label=train$is_pln,weight=train$row_weight,nthread=2),
      nrounds=120L)
    xgb.save(model,file.path(destination,"models",paste0(mode,".json")))
    probability[[mode]] <- predict(model,xgb.DMatrix(as.matrix(tune[,use,with=FALSE]),nthread=1))
  }
  grid <- rbindlist(lapply(names(probability),function(name) {
    thresholds <- if (name=="Original") .5 else c(.35,.5,.65)
    rbindlist(lapply(thresholds,function(threshold) {
      z <- tune[,.(family,is_pln,design)]
      z[, error := (probability[[name]]>=threshold) != is_pln]
      cell <- z[,.(error=mean(error)),by=.(family,design)]
      data.table(model=name,threshold=threshold,
        mean_error=mean(cell$error),max_error=max(cell$error),
        small_pln_error=mean(cell[family=="PLN" & design %in% 4:6,error]),
        complexity=match(name,c("Original","Balanced","Stratified")),
        cutoff_distance=abs(threshold-.5))
    }))
  }))
  setorder(grid,mean_error,small_pln_error,max_error,complexity,cutoff_distance)
  fwrite(grid,file.path(destination,"results","gate_validation.csv"))
  fwrite(grid[1],file.path(destination,"models","selected_gate.csv"))
  print(grid)
}

if (phase=="evaluate_new") {
  source("../reference/count_baselines.R")
  worker <- as.integer(args[2])
  workers <- as.integer(args[3])
  choice <- fread(file.path(destination,"models","selected_gate.csv"))
  old_gate <- read_frozen_model("models/pool_gate.json")
  gate <- if (choice$model=="Original") old_gate else
    read_frozen_model(file.path(destination,"models",paste0(choice$model,".json")))
  specialist_choices <- lapply(choices$family,function(family) {
    r <- fread(file.path("settings",paste0("frozen_specialist_",family,".csv")))
    list(model=read_frozen_model(file.path("models",paste0(r$model,".json"))),weight=r$weight)
  })
  adaptive_choices <- fread("settings/frozen_adaptive.csv")
  ablations_path <- "../ablations/models"
  ablation_choices <- readRDS(file.path(ablations_path,"ablation_selection_frozen.rds"))
  ablation_models <- lapply(ablation_choices$model,function(name) {
    list(model=read_frozen_model(file.path(ablations_path,paste0(name,".json"))),
         cols=readRDS(file.path(ablations_path,paste0(name,"_spec.rds")))$features)
  })

  jobs <- jobs_for("test2")
  jobs[, seed := seed+10000000L]
  guard <- as.data.table(expand.grid(design=4:6, mechanism=1:8, replicate=1:2))
  guard[, ':='(suite="NB_small",seed=851500000L+design*1000L+mechanism*10L+replicate)]
  jobs <- rbindlist(list(jobs,guard),use.names=TRUE)
  stopifnot(nrow(jobs)==240L,!anyDuplicated(jobs$seed),
            !any(jobs$seed %in% unlist(lapply(c("fit","tune","test","test2"),function(split)jobs_for(split)$seed))))
  out <- rbindlist(lapply(seq(worker,nrow(jobs),by=workers),function(i) {
    job <- jobs[i]
    sim <- if (job$suite=="NB_small")
      with(job,simulate_profile("NB",design,mechanism,seed)) else
      with(job,simulate_case(suite,design,mechanism,seed))
    f <- both_features(sim$counts,sim$n1,sim$n0)
    q <- stratified_profile(sim$counts)
    p_old <- predict(old_gate,xgb.DMatrix(matrix(q[1:23],nrow=1),nthread=1))
    x <- if (choice$model=="Stratified") q else q[1:23]
    p_new <- predict(gate,xgb.DMatrix(matrix(x,nrow=1),nthread=1))
    scores <- as.data.frame(cbind(f[,baseline_columns,drop=FALSE],
                                  count_scores(sim$counts,sim$n1,sim$n0)))
    for (index in seq_along(choices$family)) {
      family <- choices$family[index]
      morank <- family=="MoRank37"
      pair <- predict_new_old(index,f)
      branch <- specialist_choices[[index]]
      nb <- combination(predict_features(branch$model,f,features_for(family,colnames(f))),
                        pair$old,branch$weight)
      adaptive <- adaptive_choices[family==choices$family[index]]
      weight <- with(adaptive,mix_weight(p_old,w_nb,w_pln,rule))
      prefix <- if (morank) "MoRank" else "Classical"
      scores[[if (morank) "WSKR_A" else "ClassicalStack"]] <- pair$old
      scores[[choices$reported_method[index]]] <- combination(pair$new,pair$old,choices$weight[index])
      scores[[adaptive$reported_method]] <- combination(pair$new,pair$old,weight)
      scores[[paste(prefix,"final")]] <- if (p_old>=.5) pair$new else nb
      scores[[paste(prefix,"balanced routing")]] <- if (p_new>=choice$threshold) pair$new else nb
    }
    for (index in seq_len(nrow(ablation_choices))) {
      m <- ablation_models[[index]]
      scores[[ablation_choices$family[index]]] <- predict_features(m$model,f,m$cols)
    }
    hit <- function(score) {
      boundary <- sort(score,decreasing=TRUE)[35L]
      above <- score>boundary
      tied <- score==boundary
      sum(sim$truth & above)+(35-sum(above))*sum(sim$truth & tied)/sum(tied)
    }
    stopifnot(all(is.finite(as.matrix(scores))))
    cat("TEST",job$suite,job$seed,"routes",round(p_old,3),round(p_new,3),"\n")
    flush.console()
    cbind(job,data.table(n1=sim$n1,n0=sim$n0,method=names(scores),
          auc=vapply(scores,auc,0,truth=sim$truth),hits35=vapply(scores,hit,0),
          p_old=p_old,p_new=p_new,threshold=choice$threshold))
  }))
  fwrite(out,file.path(destination,"results",sprintf("evaluation_%02d.csv",worker)))
}
