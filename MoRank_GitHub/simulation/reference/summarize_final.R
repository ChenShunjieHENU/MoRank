suppressPackageStartupMessages(library(data.table));setDTthreads(1)
if(length(commandArgs(TRUE))) setwd(commandArgs(TRUE)[1])
files<-file.path('results',paste0('final_worker',1:12,'.csv'));stopifnot(all(file.exists(files)))
a<-rbindlist(lapply(files,fread));stopifnot(nrow(a)==4800,uniqueN(a,by=c('scenario','design','rep','method'))==4800,uniqueN(a$method)==20)
stopifnot(all(a[, .N, by=.(scenario,design,method)]$N==10),all(is.finite(a$auc)),all(a$auc>=0&a$auc<=1))
methods<-unique(a$method);candidate<-'WSKR_A';comparators<-setdiff(methods,candidate)
fwrite(a,'results/final_auc_all.csv')
s<-a[,.(AUC=mean(auc),SD=sd(auc),n=.N),by=.(scenario,design,n1,n0,method)];fwrite(s,'results/final_auc_by_cell.csv')
macro<-a[,{
  z<-.SD[,.(mu=mean(auc),variance=var(auc)/.N,n=.N),by=.(scenario,design)]
  mu<-mean(z$mu);se<-sqrt(sum(z$variance))/24
  df<-if(sum(z$variance)>0)sum(z$variance)^2/sum(z$variance^2/(z$n-1)) else Inf
  ci<-if(se>0)mu+c(-1,1)*qt(.975,df)*se else c(mu,mu)
  list(macro_AUC=mu,ci_low=ci[1],ci_high=ci[2],mc_se=se,n_datasets=.N)
},by=method];setorder(macro,-macro_AUC);fwrite(macro,'results/final_macro_ranking.csv')
paired<-function(z){
 n<-nrow(z);mu<-mean(z$d);se<-sd(z$d)/sqrt(n);ci<-if(se>0)mu+c(-1,1)*qt(.975,n-1)*se else c(mu,mu)
 list(difference=mu,ci_low=ci[1],ci_high=ci[2],p=if(se>0)2*pt(-abs(mu/se),n-1) else if(mu==0)1 else 0,n=n)
}
macro_rows<-cell_rows<-list()
for(mm in comparators){
 w<-dcast(a[method%in%c(candidate,mm)],scenario+design+n1+n0+rep~method,value.var='auc');w[,d:=get(candidate)-get(mm)]
 for(sc in unique(w$scenario))for(di in 1:3){z<-w[scenario==sc&design==di];r<-paired(z)
  cell_rows[[paste(sc,di,mm)]]<-data.table(scenario=sc,design=di,n1=z$n1[1],n0=z$n0[1],baseline=mm,as.data.table(r))}
 # Estimand fixes equal weights for24cells; simulation datasets are independent within each cell.
 z<-w[,.(mu=mean(d),variance=var(d)/.N,n=.N),by=.(scenario,design)]
 mu<-mean(z$mu);se<-sqrt(sum(z$variance))/24
 df<-if(sum(z$variance)>0)sum(z$variance)^2/sum(z$variance^2/(z$n-1)) else Inf
 ci<-if(se>0)mu+c(-1,1)*qt(.975,df)*se else c(mu,mu)
 macro_rows[[mm]]<-data.table(baseline=mm,difference=mu,ci_low=ci[1],ci_high=ci[2],mc_se=se,df=df,p=if(se>0)2*pt(-abs(mu/se),df) else if(mu==0)1 else 0,n_datasets=240)
}
mac<-rbindlist(macro_rows);mac[,q_BH_all19:=p.adjust(p,'BH')];fwrite(mac,'results/final_macro_paired_comparisons.csv')
cells<-rbindlist(cell_rows);cells[,q_BH_all456:=p.adjust(p,'BH')];fwrite(cells,'results/final_cell_paired_comparisons.csv')
regret<-s[,.(best_cell_AUC=max(AUC),candidate_AUC=AUC[method==candidate],best_method=method[which.max(AUC)]),by=.(scenario,design,n1,n0)]
regret[,candidate_minus_cell_best:=candidate_AUC-best_cell_AUC];fwrite(regret,'results/final_cell_regret.csv')
if(file.exists('results/stress_auc.csv')){
 st<-fread('results/stress_auc.csv');stopifnot(nrow(st)==1280,uniqueN(st,by=c('scenario','design','rep','method'))==1280)
 r<-st[,.(macro_AUC=mean(auc),n=.N),by=method];setorder(r,-macro_AUC);fwrite(r,'results/stress_macro_ranking.csv')
 fwrite(st[,.(AUC=mean(auc),SD=sd(auc),n=.N),by=.(scenario,design,n1,n0,method)],'results/stress_auc_by_cell.csv')
 ss<-list()
 for(mm in setdiff(unique(st$method),candidate)){
  w<-dcast(st[method%in%c(candidate,mm)],scenario+design+rep~method,value.var='auc');w[,d:=get(candidate)-get(mm)]
  z<-w[,.(mu=mean(d),variance=var(d)/.N,n=.N),by=.(scenario,design)];mu<-mean(z$mu);se<-sqrt(sum(z$variance))/16
  df<-sum(z$variance)^2/sum(z$variance^2/(z$n-1));ci<-mu+c(-1,1)*qt(.975,df)*se
  ss[[mm]]<-data.table(baseline=mm,difference=mu,ci_low=ci[1],ci_high=ci[2],p=2*pt(-abs(mu/se),df),n_datasets=80)
 }
 sp<-rbindlist(ss);sp[,q_BH_all15:=p.adjust(p,'BH')];fwrite(sp,'results/stress_macro_paired_comparisons.csv')
}
cat('ALL240_FINAL_DATASETS_COMPLETE;all20methods\n');print(macro);print(mac[order(difference)][1:4]);print(regret[order(candidate_minus_cell_best)][1:5])
