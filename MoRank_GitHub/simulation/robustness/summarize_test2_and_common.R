suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
expected <- c(test2=192L,neighbor_final=24L)
args <- commandArgs(TRUE)
root <- normalizePath(args[1],winslash="/"); periods <- args[-1]
if(length(periods)==0) periods <- names(expected)
for (period in periods) {
  paths <- list.files(if(period=="test2")file.path(root,"robustness") else file.path(root,"common/results"),pattern=paste0("^",period,"_[0-9]+\\.csv$"),full.names=TRUE)
  data <- rbindlist(lapply(paths,fread))
  if (period=="neighbor_final") {
    # Direct competitors are evaluated on full-precision in-memory scores.
    # Replace the AUCs derived from serialized per-gene CSV scores.
    direct_paths <- list.files(file.path(root,"common/results"),
                              "^neighbor_reference_auc_[0-9]+[.]csv$",full.names=TRUE)
    stopifnot(length(direct_paths)==6L)
    direct <- rbindlist(lapply(direct_paths,fread))
    direct[,seed:=702611060L+match(scenario,c("S1","S2","S3","S4","S6","S7","S8","S9"))*1000L+rep]
    stopifnot(!anyDuplicated(direct[,.(seed,method)]))
    data[direct,on=.(seed,method),auc:=i.auc]
  }
  stopifnot(uniqueN(data$seed)==expected[[period]],!anyDuplicated(data[,.(seed,method)]))
  reference <- if(period %in% c("test2","neighbor_final")) "MoRank final" else "MoRank optimized"
  full <- data[method==reference]
  stopifnot(nrow(full)==expected[[period]])
  # Keep suite coverage identical when comparing methods available only in some suites.
  count <- data[,.(datasets=uniqueN(seed)),by=.(suite,method)]
  coverage <- full[,.(total=.N),by=suite]
  eligible <- merge(count,coverage,by="suite")[datasets==total,.(suite,method)]
  data <- data[eligible,on=.(suite,method),nomatch=0]
  cells <- data[,.(auc=mean(auc),hits35=mean(hits35),n=.N),by=.(suite,design,mechanism,method)]
  macro <- cells[,.(macro_auc=mean(auc),top35_hits=mean(hits35),n=sum(n)),by=.(suite,method)]
  setorder(macro,suite,-macro_auc,method)
  macro[,rank:=frank(-round(macro_auc,12),ties.method="min"),by=suite]
  paired <- merge(data[method!=reference],
                  full[,.(seed,new_auc=auc,new_hits=hits35)],by="seed")
  bycell <- paired[,.(difference=mean(new_auc-auc),top35_gain=mean(new_hits-hits35),
                       variance=var(new_auc-auc)/.N,n=.N),by=.(suite,method,design,mechanism)]
  contrasts <- bycell[,{
    estimate <- mean(difference)
    se <- sqrt(sum(variance))/.N
    denom <- sum(variance^2/(n-1))
    df <- if(denom>0) sum(variance)^2/denom else Inf
    radius <- qt(.975,df)*se
    p <- if(se>0) 2*pt(-abs(estimate/se),df) else if(estimate==0) 1 else 0
    .(difference=estimate,low=estimate-radius,high=estimate+radius,p=p,
      top35_gain=mean(top35_gain),datasets=sum(n))
  },by=.(suite,method)]
  contrasts[,p_holm:=p.adjust(p,"holm")]
  fwrite(data,file.path(root,"results/final_model",paste0(period,"_all.csv")))
  fwrite(cells,file.path(root,"results/final_model",paste0(period,"_cells.csv")))
  fwrite(macro,file.path(root,"results/final_model",paste0(period,"_macro.csv")))
  if (period=="test2") {
    overall <- macro[,.(macro_auc=mean(macro_auc),top35_hits=mean(top35_hits),
                        n=sum(n),suites=.N),by=method][suites==4L]
    setorder(overall,-macro_auc,method)
    overall[,rank:=frank(-macro_auc,ties.method="min")]
    fwrite(overall,file.path(root,"results/final_model/test2_equal_suite_means.csv"))
  }
  fwrite(contrasts,file.path(root,"results/final_model",paste0(period,"_contrasts.csv")))
  print(macro[rank<=4])
  print(contrasts[method %in% c("Classical retrained","WSKR_A","ClassicalStack","waddR","MMDAgg")])
}


