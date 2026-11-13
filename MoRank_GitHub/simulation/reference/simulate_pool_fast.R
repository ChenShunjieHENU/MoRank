simulate_pool_fast <- function(scenario,n_per_group,config) {
  gene_params<-generate_parameters(config);truth<-make_truth(scenario,config)
  lib<-simulate_library(2L*n_per_group,config)
  cl<-lib[seq_len(n_per_group)];tl<-lib[n_per_group+seq_len(n_per_group)]
  counts<-matrix(0L,nrow=config$p,ncol=2L*n_per_group,
    dimnames=list(gene_params$gene,paste0('sample_',seq_len(2L*n_per_group))))
  sig<-which(truth$is_signal);nul<-which(!truth$is_signal)
  # The source assigns signal rows first; preserve its gene-major RNG order.
  stopifnot(length(sig)==0 || identical(sig,seq_along(sig)))
  for(g in sig)counts[g,]<-simulate_one_gene(truth$signal_type[g],gene_params$base_mean[g],gene_params$dispersion[g],cl,tl,config)
  if(length(nul)>0) {
    mu<-as.vector(t(outer(gene_params$base_mean[nul],c(cl,tl))))
    ph<-rep(gene_params$dispersion[nul],each=2L*n_per_group)
    v<-rnbinom(length(mu),mu=pmax(mu,1e-8),size=1/pmax(ph,1e-8))
    counts[nul,]<-matrix(v,nrow=length(nul),byrow=TRUE)
  }
  list(counts=counts,group=c(rep('control',n_per_group),rep('case',n_per_group)),truth=truth,gene_params=gene_params)
}
