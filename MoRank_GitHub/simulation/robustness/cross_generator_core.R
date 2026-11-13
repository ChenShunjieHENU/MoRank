generate_pln<-function(scenario,n1,n0,seed,G=2000L,signal_n=120L){
 stopifnot(scenario %in% 1:4,n1>1,n0>1,G>=signal_n)
 set.seed(seed);N<-n1+n0;truth<-seq_len(G)<=signal_n
 base<-pmin(pmax(rlnorm(G,log(20),1.25),.01),5000)
 sigmalog<-runif(G,.25,.75);depth<-rlnorm(N,-.35^2/2,.35);depth<-depth/mean(depth)
 latent<-matrix(rlnorm(G*N,rep(-sigmalog^2/2,N),rep(sigmalog,N)),G,N)
 ci<-n0+seq_len(n1);ni<-seq_len(n0);si<-which(truth)
 if(scenario==1){effect<-sample(c(-1,1),signal_n,replace=TRUE)*runif(signal_n,.4,1);latent[si,ci]<-latent[si,ci,drop=FALSE]*2^effect}
 if(scenario==2){
  latent[si,ni]<-latent[si,ni,drop=FALSE]*matrix(rbinom(signal_n*n0,1,.95)/.95,signal_n,n0)
  latent[si,ci]<-latent[si,ci,drop=FALSE]*matrix(rbinom(signal_n*n1,1,.55)/.55,signal_n,n1)
 }
 if(scenario==3){
  # Gamma and lognormal multipliers match each signalgene's mean and variance.
  v<-exp(sigmalog[si]^2)-1
  latent[si,ci]<-matrix(rgamma(signal_n*n1,shape=rep(1/v,n1),scale=rep(v,n1)),signal_n,n1)
 }
 if(scenario==4){
  # Three-state distributions share mean1 and variance0.1 but differ in tailmass.
  latent[si,ni]<-matrix(sample(c(.5,1,1.5),signal_n*n0,replace=TRUE,prob=c(.2,.6,.2)),signal_n,n0)
  latent[si,ci]<-matrix(sample(c(1-sqrt(.5),1,1+sqrt(.5)),signal_n*n1,replace=TRUE,prob=c(.1,.8,.1)),signal_n,n1)
 }
 mu<-latent*base;mu<-sweep(mu,2,depth,'*')
 counts<-matrix(rpois(G*N,pmax(mu,1e-8)),G,N,dimnames=list(paste0('gene_',seq_len(G)),paste0('sample_',seq_len(N))))
 list(counts=counts,truth=truth,n1=n1,n0=n0,seed=seed,scenario=scenario)
}
