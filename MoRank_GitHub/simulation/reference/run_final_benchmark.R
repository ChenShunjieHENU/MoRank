suppressPackageStartupMessages(library(data.table));setDTthreads(1)
suppressPackageStartupMessages(library(xgboost))
source('original_simulation.R');source('simulate_pool_fast.R');source('redesign_core.R');source('count_baselines.R')
source('model_io.R')
stopifnot(file.exists('FINAL_MODEL_FREEZE.json'))
args<-commandArgs(TRUE);worker<-as.integer(args[1]);workers<-if(length(args)>=2L)as.integer(args[2]) else 12L;stopifnot(worker%in%1:workers)
output <- if(length(args)>=3L) args[3] else '.'
dir.create(file.path(output,'final_raw'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(output,'results'),recursive=TRUE,showWarnings=FALSE)
sel<-readRDS('models/frozen_selection.rds')$selected
learners<-lapply(seq_len(nrow(sel)),function(i){nm<-sel$model[i];list(family=sel$family[i],model=read_frozen_model(file.path('models',paste0(nm,'.json'))),spec=readRDS(file.path('models',paste0(nm,'_spec.rds'))))})
scenarios<-c('S1','S2','S3','S4','S6','S7','S8','S9');designs<-list(c(50L,50L),c(200L,200L),c(541L,59L));cfg<-default_config();start<-Sys.time();metrics<-list()
for(di in 1:3)for(si in 1:8){
 cell<-(di-1)*8+si;if((cell-1)%%workers+1!=worker)next
 nn<-designs[[di]];pooln<-max(nn)
 for(ri in seq_len(if(length(args)>=4L)as.integer(args[4]) else 10L)){
   seed<-202611060L+3L*1000000L+di*100000L+si*1000L+ri;set.seed(seed)
   sim<-simulate_pool_fast(scenarios[si],pooln,cfg);ix<-c(seq_len(nn[2]),pooln+seq_len(nn[1]));counts<-sim$counts[,ix,drop=FALSE]
   lib<-colSums(counts);lib[lib<=0]<-1;x<-log2(t(t(counts)/lib)*1e6+1)
   f<-distribution_features(x,nn[2]+seq_len(nn[1]),seq_len(nn[2]));scores<-f[,baseline_columns,drop=FALSE]
   for(learner in learners)scores<-cbind(scores,setNames(data.frame(predict(learner$model,xgb.DMatrix(f[,learner$spec$features,drop=FALSE],nthread=1))),learner$family))
   scores<-cbind(scores,count_scores(counts,nn[1],nn[2]));scores<-as.matrix(scores)
   stopifnot(nrow(scores)==2000,all(is.finite(scores)))
   truth<-sim$truth$is_signal;au<-vapply(seq_len(ncol(scores)),function(i)auc(scores[,i],truth),0)
   metrics[[length(metrics)+1]]<-data.table(cell=cell,scenario=scenarios[si],design=di,n1=nn[1],n0=nn[2],rep=ri,seed=seed,method=colnames(scores),auc=au)
   saveRDS(list(scenario=scenarios[si],design=di,n1=nn[1],n0=nn[2],rep=ri,seed=seed,features=f,scores=scores,truth=truth,
    counts=if(ri==1)counts else NULL),file.path(output,'final_raw',sprintf('d%d_s%d_r%02d.rds',di,si,ri)))
   fwrite(rbindlist(metrics),file.path(output,'results',paste0('final_worker',worker,'.csv')))
   cat('worker',worker,'cell',cell,'rep',ri,'COMPLETE elapsed',round(as.numeric(difftime(Sys.time(),start,units='secs')),1),'\n');flush.console()
 }
}
writeLines(capture.output(sessionInfo()),file.path(output,'results',paste0('final_worker',worker,'_R_sessionInfo.txt')))
cat('FINAL_WORKER_COMPLETE',worker,'datasets',length(metrics),'\n')
