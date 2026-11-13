suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
R <- normalizePath(commandArgs(TRUE)[1], winslash="/")
summary_matrix<-function(a,prefix,adjust='holm'){
 stopifnot(!anyDuplicated(a[,.(scenario,design,rep,method)]),all(is.finite(a$auc)))
 ncell<-uniqueN(a,by=c('scenario','design'));candidate<-'WSKR_A'
 rank<-a[,{
  z<-.SD[,.(mu=mean(auc),v=var(auc)/.N,n=.N),by=.(scenario,design)]
  se<-sqrt(sum(z$v))/ncell;df<-sum(z$v)^2/sum(z$v^2/(z$n-1));mu<-mean(z$mu)
  list(macro_AUC=mu,ci_low=mu-qt(.975,df)*se,ci_high=mu+qt(.975,df)*se,mc_se=se,datasets=.N)
 },by=method];setorder(rank,-macro_AUC);fwrite(rank,file.path(R,'results',paste0(prefix,'_macro.csv')))
 contrasts<-list()
 for(mm in setdiff(unique(a$method),candidate)){
  w<-dcast(a[method %in% c(candidate,mm)],scenario+design+rep~method,value.var='auc');stopifnot(!anyNA(w));w[,d:=get(candidate)-get(mm)]
  z<-w[,.(mu=mean(d),v=var(d)/.N,n=.N),by=.(scenario,design)]
  se<-sqrt(sum(z$v))/ncell;df<-if(sum(z$v)>0)sum(z$v)^2/sum(z$v^2/(z$n-1)) else Inf;mu<-mean(z$mu)
  contrasts[[mm]]<-data.table(baseline=mm,difference=mu,ci_low=mu-qt(.975,df)*se,ci_high=mu+qt(.975,df)*se,
    p=if(se>0)2*pt(-abs(mu/se),df) else if(mu==0)1 else 0,n_datasets=nrow(w))
 }
 cc<-rbindlist(contrasts);cc[,q_Holm:=p.adjust(p,adjust)];fwrite(cc,file.path(R,'results',paste0(prefix,'_contrasts.csv')))
 fwrite(a[,.(AUC=mean(auc),SD=sd(auc),n=.N),by=.(scenario,design,method)],file.path(R,'results',paste0(prefix,'_by_cell.csv')))
 print(rank);print(cc)
}
a <- fread(file.path(R, "results/final_auc_all.csv"))[method=="WSKR_A"]
b <- fread(file.path(R, "results/ablation_oldfinal_auc.csv"))
stopifnot(nrow(a)==240L,nrow(b)==960L)
summary_matrix(rbindlist(list(a,b),use.names=TRUE,fill=TRUE),"ablation")
