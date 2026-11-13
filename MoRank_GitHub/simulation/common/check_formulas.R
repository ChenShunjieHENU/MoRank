suppressPackageStartupMessages({library(data.table);library(jsonlite)})
R<-normalizePath(commandArgs(TRUE)[1],winslash='/')
# Exhaust all15 label allocations, independent of the model's training code.
sets<-list(tied=c(0,0,1,2,2,8),continuous=c(-2,-.2,.4,1,2,4));out<-list();ix<-0L
for(pool in names(sets)){
 x<-sets[[pool]];N<-length(x);n1<-2;n0<-N-n1;ne<-n1*n0/N;H<-(rank(x,ties.method='average')-.5)/N
 fs<-list();pr<-rep(1,N);cu<-2*H-1
 for(k in 1:8){if(k>1){nx<-((2*k-1)*(2*H-1)*cu-(k-1)*pr)/k;pr<-cu;cu<-nx};fs[[paste0('P',k)]]<-cu}
 for(q in c(.05,.1,.2,.3,.5,.7,.8,.9,.95))fs[[paste0('CDF',q)]]<-as.numeric(H<=q)
 for(fn in names(fs)){
  y<-fs[[fn]];comb<-combn(N,n1);delta<-apply(comb,2,function(ii)mean(y[ii])-mean(y[-ii]));sf<-var(y)
  qf<-pmin(1e6,ne*delta^2/max(sf,1e-12));err<-abs(mean(delta^2)-sf/ne)
  stopifnot(abs(mean(delta))<1e-12,err<1e-12,mean(qf)<=1+1e-10,max(qf)<=N-1+1e-10)
  ix<-ix+1L;out[[ix]]<-data.table(pool=pool,basis=fn,N=N,n1=n1,allocations=ncol(comb),mean_delta=mean(delta),var_exact=mean(delta^2),var_formula=sf/ne,max_abs_error=err,mean_Q=mean(qf),max_Q=max(qf),status='PASS')
 }
}
fwrite(rbindlist(out),file.path(R,'results/formula_rank_enumeration.csv'),bom=TRUE)
mom<-function(x,w){mu<-sum(w*x);v<-sum(w*(x-mu)^2);z<-(x-mu)/sqrt(v);S<-sum(w*z^3);K<-sum(w*z^4)-3;c(S=S,K=K,T=sign(K)*log1p(abs(K)))}
x<-c(-2,-1,.2,.5,1,1.1,2,7);w<-rep(1/length(x),length(x));m<-mom(x,w);z<-(x-mean(x))/sqrt(mean((x-mean(x))^2));S<-m['S'];K<-m['K'];psiS<-z^3-3*z-S-1.5*S*(z^2-1);psiK<-z^4-4*S*z-(K+3)-2*(K+3)*(z^2-1);ps<-cbind(S=psiS,K=psiK,T=psiK/(1+abs(K)));eps<-1e-7
nf<-t(vapply(seq_along(x),function(i){ww<-(1-eps)*w;ww[i]<-ww[i]+eps;(mom(x,ww)-m)/eps},numeric(3)));er<-abs(nf-ps);stopifnot(max(er)<1e-4)
zz<-rbindlist(lapply(1:3,function(k)data.table(point=seq_along(x),component=colnames(ps)[k],formula=ps[,k],finite_difference=nf[,k],abs_error=er[,k],epsilon=eps,status='PASS')));fwrite(zz,file.path(R,'results/formula_influence_finite_difference.csv'),bom=TRUE)
write_json(list(status='PASS',scope='Finite exact-randomization examples and numerical influence derivatives; no global optimality proof.',rank_cases=length(out),allocations_per_case=15,influence_derivatives=nrow(zz),max_influence_error=max(er),equations_changed=FALSE),file.path(R,'qa/formula_verification.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
cat('FORMULA_CHECK_PASS',length(out),'rank cases;',nrow(zz),'influence derivatives\n')


