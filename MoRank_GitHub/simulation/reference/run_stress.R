suppressPackageStartupMessages(library(data.table));setDTthreads(1);suppressPackageStartupMessages(library(xgboost))
source('original_simulation.R');source('simulate_pool_fast.R');source('redesign_core.R');stopifnot(file.exists('FINAL_MODEL_FREEZE.json'))
source('model_io.R')
args <- commandArgs(TRUE)
output <- if(length(args)) args[1] else '.'
dir.create(file.path(output,'stress_raw'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(output,'results'),recursive=TRUE,showWarnings=FALSE)
sel<-readRDS('models/frozen_selection.rds')$selected
learners<-lapply(seq_len(nrow(sel)),function(i){nm<-sel$model[i];list(family=sel$family[i],model=read_frozen_model(file.path('models',paste0(nm,'.json'))),spec=readRDS(file.path('models',paste0(nm,'_spec.rds'))))})
scenarios<-c('S1','S2','S3','S4','S6','S7','S8','S9');designs<-list(c(50L,50L),c(541L,59L));cfg<-default_config()
cfg$tail_pi<-.02;cfg$tail_ratio<-8;cfg$dispersion_ratio<-1.4;cfg$skew_sigma_case<-.4;cfg$proportion_case<-.5
cfg$modality_low_ratio<-.75;cfg$modality_high_ratio<-1.25;cfg$mean_log2fc_min<-.3;cfg$mean_log2fc_max<-.8
metrics<-list();st<-Sys.time()
for(di in 1:2)for(si in 1:8)for(ri in 1:5){
 seed<-202611060L+4L*1000000L+di*100000L+si*1000L+ri;set.seed(seed);nn<-designs[[di]];pooln<-max(nn)
 sim<-simulate_pool_fast(scenarios[si],pooln,cfg);ix<-c(seq_len(nn[2]),pooln+seq_len(nn[1]));counts<-sim$counts[,ix,drop=FALSE]
 lib<-colSums(counts);lib[lib<=0]<-1;x<-log2(t(t(counts)/lib)*1e6+1);f<-distribution_features(x,nn[2]+seq_len(nn[1]),seq_len(nn[2]));scores<-f[,baseline_columns,drop=FALSE]
 for(learner in learners)scores<-cbind(scores,setNames(data.frame(predict(learner$model,xgb.DMatrix(f[,learner$spec$features,drop=FALSE],nthread=1))),learner$family))
 scores<-as.matrix(scores);truth<-sim$truth$is_signal
 metrics[[length(metrics)+1]]<-data.table(scenario=scenarios[si],design=di,n1=nn[1],n0=nn[2],rep=ri,seed=seed,method=colnames(scores),auc=vapply(seq_len(ncol(scores)),function(i)auc(scores[,i],truth),0))
 saveRDS(list(seed=seed,scores=scores,truth=truth),file.path(output,'stress_raw',sprintf('d%d_s%d_r%02d.rds',di,si,ri)))
 if(ri==5){fwrite(rbindlist(metrics),file.path(output,'results/stress_auc.csv'));cat('STRESS design',di,'scenario',scenarios[si],'COMPLETE elapsed',round(as.numeric(difftime(Sys.time(),st,units='secs')),1),'\n');flush.console()}
}
cat('STRESS_ALL80_COMPLETE\n')
