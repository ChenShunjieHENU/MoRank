suppressPackageStartupMessages({library(data.table);library(xgboost)})
setDTthreads(1)
old <- normalizePath('../reference')
R <- normalizePath(Sys.getenv('MORANK_TRAIN_OUTPUT','refit'),mustWork=FALSE)
dir.create(file.path(R,'models'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(R,'results'),recursive=TRUE,showWarnings=FALSE)
development <- Sys.getenv('MORANK_DEVELOPMENT',file.path(old,'cache'))
setwd(old);source('redesign_core.R');source('model_io.R')
if(!file.exists(file.path(R,'models/ablation_selection_frozen.rds'))){
fit<-lapply(sort(list.files(file.path(development,'fit'),full.names=TRUE)),readRDS)
val<-lapply(sort(list.files(file.path(development,'validation'),full.names=TRUE)),readRDS)
stopifnot(length(fit)==144,length(val)==72)
x<-do.call(rbind,lapply(fit,`[[`,'features'));y<-as.integer(unlist(lapply(fit,`[[`,'truth')))
set.seed(202611066);o<-sample.int(length(y));x<-x[o,];y<-y[o]
allcols<-colnames(x);ranks<-grep('^(Rank|CDF)',allcols,value=TRUE)
moments<-c('OriginalWSKR','SR','KR','Decancel','ShapeNoise','SkewNoise','TailNoise')
variants<-list(DropOriginalWSKR=setdiff(allcols,'OriginalWSKR'),DropMomentBlock=setdiff(allcols,moments),
 DropRankBlock=setdiff(allcols,ranks),RankOnly=c(ranks,'log_n_eff','imbalance','tie_fraction'))
stopifnot(lengths(variants)==c(36,30,17,23))
grid<-expand.grid(max_depth=c(2L,3L,4L),nrounds=c(150L,400L));records<-list();selected<-list();started<-Sys.time()
for(fam in commandArgs(TRUE)[1]){
 cols<-variants[[fam]];dm<-xgb.DMatrix(x[,cols,drop=FALSE],label=y,nthread=1)
 for(i in seq_len(nrow(grid))){
  params<-list(objective='binary:logistic',eval_metric='logloss',tree_method='hist',nthread=1,
   max_depth=grid$max_depth[i],eta=.05,min_child_weight=20,lambda=10,subsample=.8,colsample_bytree=.9,seed=202611067)
  mod<-xgb.train(params=params,data=dm,nrounds=grid$nrounds[i],verbose=0)
  va<-vapply(val,function(z)auc(predict(mod,xgb.DMatrix(z$features[,cols,drop=FALSE],nthread=1)),z$truth),0.)
  nm<-paste0(fam,'_d',grid$max_depth[i],'_t',grid$nrounds[i])
  records[[nm]]<-data.table(family=fam,model=nm,features=length(cols),max_depth=grid$max_depth[i],nrounds=grid$nrounds[i],validation_macro_AUC=mean(va))
  xgb.save(mod,file.path(R,'models',paste0(nm,'.json')))
  saveRDS(list(features=cols,params=params,nrounds=grid$nrounds[i]),file.path(R,'models',paste0(nm,'_spec.rds')))
 }
 fam_name<-fam;sel<-rbindlist(records)[family==fam_name];setorder(sel,-validation_macro_AUC,max_depth,nrounds);selected[[fam]]<-sel[1]
 cat(fam,'SELECTED',sel$model[1],'features',length(cols),'elapsed',round(as.numeric(difftime(Sys.time(),started,units='secs'))),'\n');flush.console()
}
selection<-rbindlist(selected);fwrite(rbindlist(records),file.path(R,'results',paste0('grid_',fam,'.csv')));fwrite(selection,file.path(R,'results',paste0('selected_',fam,'.csv')))
}
