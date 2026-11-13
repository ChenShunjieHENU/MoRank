suppressPackageStartupMessages({library(data.table);library(xgboost)})
setDTthreads(1)
R <- normalizePath('.')
old <- normalizePath('../reference')
setwd(old);source('redesign_core.R');source('model_io.R')
selection<-readRDS(file.path(R,'models/ablation_selection_frozen.rds'))
learners<-lapply(seq_len(nrow(selection)),function(i){nm<-selection$model[i];list(family=selection$family[i],model=read_frozen_model(file.path(R,'models',paste0(nm,'.json'))),spec=readRDS(file.path(R,'models',paste0(nm,'_spec.rds'))))})
dir.create(Sys.getenv('MORANK_AUDIT_RESULTS'),recursive=TRUE,showWarnings=FALSE)
metrics<-list();preds<-list();files<-sort(list.files(file.path(Sys.getenv('MORANK_PRIMARY_DIR'),'final_raw'),pattern='[.]rds$',full.names=TRUE));stopifnot(length(files)==240)
for(j in seq_along(files)){
 z<-readRDS(files[j]);f<-z$features
 for(a in learners){
  sc<-predict(a$model,xgb.DMatrix(f[,a$spec$features,drop=FALSE],nthread=1))
  metrics[[length(metrics)+1]]<-data.table(scenario=z$scenario,design=z$design,n1=z$n1,n0=z$n0,rep=z$rep,seed=z$seed,method=a$family,auc=auc(sc,z$truth))
  preds[[length(preds)+1]]<-data.table(scenario=z$scenario,design=z$design,rep=z$rep,gene_index=seq_along(sc),truth=z$truth,method=a$family,score=sc)
 }
 if(j%%40==0){cat('Evaluated oldfinal',j,'/240\n');flush.console()}
}
fwrite(rbindlist(metrics),file.path(Sys.getenv('MORANK_AUDIT_RESULTS'),'ablation_oldfinal_auc.csv'))
fwrite(rbindlist(preds),file.path(Sys.getenv('MORANK_AUDIT_RESULTS'),'ablation_oldfinal_predictions.csv.gz'))
cat('ABLATIONS_COMPLETED\n')
