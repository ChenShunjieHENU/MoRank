R <- normalizePath(Sys.getenv('MORANK_COMMON_DIR'),mustWork=TRUE)
old <- normalizePath('../reference')
suppressPackageStartupMessages({library(data.table);library(xgboost)});setDTthreads(1)
setwd(old);source('original_simulation.R');source('simulate_pool_fast.R');source('redesign_core.R');source('model_io.R');source('count_baselines.R')
dir.create(file.path(R,'data/neighbor'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(R,'results'),showWarnings=FALSE)
sel<-readRDS(file.path(old,'models/frozen_selection.rds'))$selected
learners<-lapply(seq_len(nrow(sel)),function(i){nm<-sel$model[i];list(family=sel$family[i],model=read_frozen_model(file.path(old,'models',paste0(nm,'.json'))),spec=readRDS(file.path(old,'models',paste0(nm,'_spec.rds'))))})
scenarios<-c('S1','S2','S3','S4','S6','S7','S8','S9');cfg<-default_config();metrics<-list();start<-Sys.time()
for(si in 1:8)for(ri in 1:3){
 id<-(si-1)*3+ri;if((id-1)%%6+1!=as.integer(commandArgs(TRUE)[1]))next;seed<-702611060L+si*1000L+ri;set.seed(seed);sim<-simulate_pool_fast(scenarios[si],50L,cfg)
 x<-log2(sweep(sim$counts,2,colSums(sim$counts),'/')*1e6+1);f<-distribution_features(x,51:100,1:50)
 score<-f[,baseline_columns,drop=FALSE]
 for(a in learners)score<-cbind(score,setNames(data.frame(predict(a$model,xgb.DMatrix(f[,a$spec$features,drop=FALSE],nthread=1))),a$family))
 score<-as.matrix(cbind(score,count_scores(sim$counts,50,50)))
 set.seed(seed+10000000L);ix<-sample(c(sample(which(sim$truth$is_signal),80),sample(which(!sim$truth$is_signal),320)))
 out<-data.table(dataset=id,scenario=scenarios[si],design=1L,n1=50L,n0=50L,rep=ri,seed=seed,gene_index=ix,truth=sim$truth$is_signal[ix])
 out<-cbind(out,as.data.table(score[ix,,drop=FALSE]));fwrite(out,file.path(R,'data/neighbor',sprintf('scores_%02d.csv',id)))
 fwrite(as.data.table(x[ix,,drop=FALSE]),file.path(R,'data/neighbor',sprintf('expression_%02d.csv.gz',id)))
 for(mm in colnames(score))metrics[[length(metrics)+1]]<-data.table(dataset=id,scenario=scenarios[si],design=1L,rep=ri,method=mm,auc=auc(score[ix,mm],sim$truth$is_signal[ix]))
 fwrite(rbindlist(metrics),file.path(R,'results',sprintf('neighbor_reference_auc_%02d.csv',as.integer(commandArgs(TRUE)[1]))))
 cat('NEIGHBOR_DATA_COMPLETE',id,'/24 elapsed',round(as.numeric(difftime(Sys.time(),start,units='secs'))),'\n');flush.console()
}
cat('ALL_NEIGHBOR_DATA_COMPLETE\n')

