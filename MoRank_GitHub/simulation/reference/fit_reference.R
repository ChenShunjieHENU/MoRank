suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(xgboost));source('redesign_core.R')
output <- Sys.getenv('MORANK_TRAIN_OUTPUT', 'refit')
dir.create(file.path(output,'models'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(output,'results'),recursive=TRUE,showWarnings=FALSE)
fitfiles<-sort(list.files(file.path(Sys.getenv('MORANK_DEVELOPMENT','cache'),'fit'),full.names=TRUE));valfiles<-sort(list.files(file.path(Sys.getenv('MORANK_DEVELOPMENT','cache'),'validation'),full.names=TRUE))
stopifnot(length(fitfiles)==144,length(valfiles)==72)
fit<-lapply(fitfiles,readRDS);val<-lapply(valfiles,readRDS)
x<-do.call(rbind,lapply(fit,`[[`,'features'));y<-as.integer(unlist(lapply(fit,`[[`,'truth')))
stopifnot(nrow(x)==length(y),!anyNA(x));set.seed(202611066);o<-sample.int(length(y));x<-x[o,];y<-y[o]
timing<-Sys.time();logs<-list();allauc<-list();best<-list();grid<-expand.grid(max_depth=c(2L,3L,4L),nrounds=c(150L,400L))
eval_scores<-function(scorefun,nm){rbindlist(lapply(val,function(u)data.table(method=nm,scenario=u$scenario,design=u$design,rep=u$rep,seed=u$seed,auc=auc(scorefun(u),u$truth))))}
for(nm in baseline_columns)allauc[[nm]]<-eval_scores(function(u)u$features[,nm],nm)
for(family in commandArgs(TRUE)[1]){
 cols<-if(family=='ClassicalStack')classic_features else colnames(x)
 dm<-xgb.DMatrix(x[,cols,drop=FALSE],label=y,nthread=1)
 for(i in seq_len(nrow(grid))){
   par<-list(objective='binary:logistic',eval_metric='logloss',tree_method='hist',nthread=1,
    max_depth=grid$max_depth[i],eta=.05,min_child_weight=20,lambda=10,subsample=.8,colsample_bytree=.9,seed=202611067)
   model<-xgb.train(params=par,data=dm,nrounds=grid$nrounds[i],verbose=0)
   nm<-paste0(family,'_d',grid$max_depth[i],'_t',grid$nrounds[i]);scorefun<-function(u)predict(model,xgb.DMatrix(u$features[,cols,drop=FALSE],nthread=1))
   a<-eval_scores(scorefun,nm);mac<-mean(a$auc);logs[[nm]]<-data.table(family=family,model=nm,max_depth=grid$max_depth[i],nrounds=grid$nrounds[i],validation_macro_AUC=mac)
   allauc[[nm]]<-a;xgb.save(model,file.path(output,'models',paste0(nm,'.json')))
   saveRDS(list(features=cols,params=par,nrounds=grid$nrounds[i],model_name=nm),file.path(output,'models',paste0(nm,'_spec.rds')))
   cat(nm,'validation_macro_AUC',format(mac,digits=8),'seconds',round(as.numeric(difftime(Sys.time(),timing,units='secs')),1),'\n');flush.console()
 }
 family_name<-family;z<-rbindlist(logs)[family==family_name];setorder(z,-validation_macro_AUC,max_depth,nrounds);best[[family]]<-z[1]
}
log<-rbindlist(logs);fwrite(log,file.path(output,'results',paste0('grid_',family,'.csv')))
v<-rbindlist(allauc);fwrite(v,file.path(output,'results',paste0('validation_',family,'.csv')))
selection<-rbindlist(best);fwrite(selection,file.path(output,'results',paste0('selection_',family,'.csv')))
baseline<-v[,.(macro_AUC=mean(auc)),by=method];setorder(baseline,-macro_AUC);fwrite(baseline,file.path(output,'results',paste0('ranking_',family,'.csv')))
saveRDS(list(selected=selection,feature_names=colnames(x),train_seeds=vapply(fit,`[[`,0,'seed'),validation_seeds=vapply(val,`[[`,0,'seed'),never_used='Final split_id3 and stress split_id4'), file.path(output,'models',paste0('selection_',family,'.rds')))
writeLines(capture.output(sessionInfo()),file.path(output,'results','R_model_sessionInfo.txt'))
cat('MODEL_SELECTION_FROZEN','\n');print(selection)

