source('original_simulation.R');source('simulate_pool_fast.R');source('redesign_core.R')
args<-commandArgs(TRUE);split<-args[1];sid<-switch(split,fit=1L,validation=2L);reps<-if(split=='fit')6L else 3L
stopifnot(!is.null(sid));dir.create(file.path(Sys.getenv('MORANK_DEVELOPMENT', 'cache'),split),recursive=TRUE,showWarnings=FALSE)
scenarios<-c('S1','S2','S3','S4','S6','S7','S8','S9');designs<-list(c(50L,50L),c(200L,200L),c(541L,59L));cfg<-default_config();start<-Sys.time()
for(di in 1:3)for(si in 1:8)for(ri in seq_len(reps)){
 
 seed<-202611060L+sid*1000000L+di*100000L+si*1000L+ri
 nn<-designs[[di]];pooln<-max(nn);set.seed(seed);sim<-simulate_pool_fast(scenarios[si],pooln,cfg)
 ix<-c(seq_len(nn[2]),pooln+seq_len(nn[1]));counts<-sim$counts[,ix,drop=FALSE]
 lib<-colSums(counts);lib[lib<=0]<-1;x<-log2(t(t(counts)/lib)*1e6+1)
 features<-distribution_features(x,nn[2]+seq_len(nn[1]),seq_len(nn[2]))
 # Feature matrices exclude truth, gene index/name, scenario and all generator parameters.
 unit<-list(split=split,split_id=sid,scenario=scenarios[si],design=di,n1=nn[1],n0=nn[2],rep=ri,seed=seed,
  features=features,truth=sim$truth$is_signal)
 saveRDS(unit,file.path(Sys.getenv('MORANK_DEVELOPMENT', 'cache'),split,sprintf('d%d_s%d_r%02d.rds',di,si,ri)))
 if(ri==reps){cat(split,'design',di,'scenario',scenarios[si],'completed',ri,'elapsed',round(as.numeric(difftime(Sys.time(),start,units='secs')),1),'\n');flush.console()}
}
cat('DEVELOPMENT_SPLIT_COMPLETE',split,'seconds',as.numeric(difftime(Sys.time(),start,units='secs')),'\n')
