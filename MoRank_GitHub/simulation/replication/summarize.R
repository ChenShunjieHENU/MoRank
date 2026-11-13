suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
paths <- list.files("results",pattern="^evaluation_[0-9]+\\.csv$",full.names=TRUE)
stopifnot(length(paths)==20L)
records <- rbindlist(lapply(paths,fread))
suites <- c("NB","Sensitivity","PLN","PLN_small")
units <- unique(records[,.(suite,design,mechanism,replicate,seed,n1,n0,p_old,p_new,threshold)])
stopifnot(nrow(units)==240L,!anyDuplicated(units$seed),
          nrow(records)==240L*32L,!anyDuplicated(records[,.(seed,method)]),
          all(is.finite(records$auc)),all(records$auc>=0 & records$auc<=1),
          all(records$hits35>=0 & records$hits35<=35))
stopifnot(identical(as.integer(units[,.(N=.N),by=suite][match(c(suites,"NB_small"),suite),N]),
                    c(72L,48L,36L,36L,48L)))

rankable <- records
selected <- fread("models/selected_gate.csv")
if (selected$model=="Original" && selected$threshold==.5) {
  for (family in c("MoRank","Classical")) {
    base <- records[method==paste(family,"final"),.(seed,auc,hits35)][order(seed)]
    alias <- records[method==paste(family,"balanced routing"),.(seed,auc,hits35)][order(seed)]
    stopifnot(identical(base,alias))
  }
  rankable <- records[!method %in% c("MoRank balanced routing","Classical balanced routing")]
}
cells <- rankable[,.(auc=mean(auc),hits35=mean(hits35),n=.N),by=.(suite,design,mechanism,method)]
macro <- cells[,.(macro_auc=mean(auc),top35_hits=mean(hits35),n=sum(n)),by=.(suite,method)]
macro[,rank:=frank(-macro_auc,ties.method="min"),by=suite]
overall <- macro[suite %in% suites,.(macro_auc=mean(macro_auc),top35_hits=mean(top35_hits),n=sum(n)),by=method]
overall[,rank:=frank(-macro_auc,ties.method="min")]
setorder(macro,suite,rank,method)
setorder(overall,rank,method)
fwrite(records,"results/independent_scores.csv")
fwrite(cells,"results/independent_cells.csv")
fwrite(macro,"results/independent_suite_means.csv")
fwrite(overall,"results/independent_equal_suite_means.csv")

units[, true_pln := startsWith(suite,"PLN")]
routes <- units[,.(datasets=.N,previous_wrong=sum((p_old>=.5)!=true_pln),
                   revised_wrong=sum((p_new>=threshold)!=true_pln)),by=.(suite,design,n1,n0)]
setorder(routes,suite,design)
fwrite(routes,"results/routing_accuracy.csv")
fwrite(units,"results/independent_seed_manifest.csv")

targets <- c("MoRank final","MoRank balanced routing")
comparators <- c("Classical final","Classical balanced routing","MoRank optimized","SR","DESeq2_BH","edgeR_BH","WSKR_A")
grid <- CJ(target=intersect(targets,rankable$method), reference=intersect(comparators,rankable$method))
paired <- rbindlist(lapply(seq_len(nrow(grid)),function(i) {
  label <- grid[i]
  first <- records[method==label$target,.(suite,design,mechanism,seed,auc, hits35)]
  second <- records[method==label$reference,.(seed,auc, hits35)]
  by_dataset <- merge(first,second,by="seed",suffixes=c("_model","_reference"))
  stopifnot(nrow(by_dataset)==240L)
  by_dataset[,delta:=auc_model-auc_reference]
  cells <- by_dataset[,.(d=mean(delta),v=var(delta)/.N,df=.N-1L,
                         hits_delta=mean(hits35_model-hits35_reference)),
                       by=.(suite,design,mechanism)]
  estimates <- cells[, {
    n <- .N
    d <- mean(d)
    se <- sqrt(sum(v))/n
    df_eff <- if (sum(v)>0) sum(v)^2/sum(v^2/df) else Inf
    half_width <- qt(.975,df_eff)*se
    .(auc_difference=d,ci_low=d-half_width,ci_high=d+half_width,
      p=if(se>0) 2*pt(-abs(d/se),df_eff) else as.numeric(d==0),top35_gain=mean(hits_delta))
  },by=suite]
  cbind(label[rep(1L,nrow(estimates))],estimates)
}))
paired[,holm_p:=p.adjust(p,"holm"),by=.(suite,target)]
fwrite(paired,"results/independent_paired_contrasts.csv")
print(overall[rank<=7])
print(macro[method %in% c(targets,"MoRank optimized","Classical final")])
print(routes)
