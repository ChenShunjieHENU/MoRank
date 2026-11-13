suppressPackageStartupMessages(library(matrixStats))
source('analysis_core.R')
positive_ratio<-function(num,den){z<-num/pmax(den,1e-12);pmin(pmax(z,0),1e6)}
moment_noise<-function(x){
 n<-ncol(x);mu<-rowMeans(x);z<-sweep(x,1,mu,'-');s<-sqrt(rowMeans(z*z));bad<-s<=1e-8;s[bad]<-1;z<-sweep(z,1,s,'/');z2<-z*z
 S<-rowMeans(z*z2);K<-rowMeans(z2*z2)-3;S[bad]<-0;K[bad]<-0;T<-sign(K)*log1p(abs(K))
 iS<-z*z2-3*z-S-1.5*S*(z2-1)
 iK<-z2*z2-4*S*z-(K+3)-2*(K+3)*(z2-1)
 iT<-iK/(1+abs(K));iS[bad,]<-0;iT[bad,]<-0
 list(S=S,K=K,T=T,vS=matrixStats::rowVars(iS)/n,vT=matrixStats::rowVars(iT)/n)
}
distribution_features<-function(x,case,control){
 x<-as.matrix(x);G<-nrow(x);N<-ncol(x);n1<-length(case);n0<-length(control);ne<-n1*n0/N
 stopifnot(N==n1+n0,n1>1,n0>1,length(intersect(case,control))==0,all(is.finite(x)))
 a<-x[,case,drop=FALSE];b<-x[,control,drop=FALSE]
 rp<-matrixStats::rowRanks(x,ties.method='average',useNames=FALSE)
 rx<-matrixStats::rowRanks(x,ties.method='max',useNames=FALSE)
 r1<-matrixStats::rowRanks(a,ties.method='average',useNames=FALSE)-.5
 r0<-matrixStats::rowRanks(b,ties.method='average',useNames=FALSE)-.5
 # Mid-CDF differences at every pooled point; repeated ties carry their multiplicity.
 d<-matrix(0,G,N);d[,case]<-r1/n1-(rp[,case,drop=FALSE]-.5-r1)/n0
 d[,control]<-(rp[,control,drop=FALSE]-.5-r0)/n1-r0/n0
 H<-(rp-.5)/N;ties<-2*(rx-rp)+1
 den<-H*(1-H)-ties/(4*N);adterms<-d*d/pmax(den,1e-12);adterms[den<=0]<-0
 ad<-(N-1)/N*ne*rowMeans(adterms);cvm<-ne*rowMeans(d*d)
 # Right-CDF differences give the exact empirical KS distance, including ties.
 rr1<-matrixStats::rowRanks(a,ties.method='max',useNames=FALSE)
 rr0<-matrixStats::rowRanks(b,ties.method='max',useNames=FALSE)
 dr<-matrix(0,G,N);dr[,case]<-rr1/n1-(rx[,case,drop=FALSE]-rr1)/n0
 dr[,control]<-(rx[,control,drop=FALSE]-rr0)/n1-rr0/n0
 ks<-matrixStats::rowMaxs(abs(dr))
 wp<-wilcox_rows(x,case,control)
 pooled_var<-matrixStats::rowVars(x);welch<-positive_ratio((rowMeans(a)-rowMeans(b))^2,matrixStats::rowVars(a)/n1+matrixStats::rowVars(b)/n0)
 prob<-seq(.01,.99,length.out=99)
 qa<-matrix(matrixStats::rowQuantiles(a,probs=prob,type=8,useNames=FALSE),nrow=G);qb<-matrix(matrixStats::rowQuantiles(b,probs=prob,type=8,useNames=FALSE),nrow=G)
 wass<-positive_ratio(ne*rowMeans((qa-qb)^2),pooled_var)
 qac<-qa-rowMeans(qa);qbc<-qb-rowMeans(qb)
 rho<-rowSums(qac*qbc)/sqrt(rowSums(qac*qac)*rowSums(qbc*qbc));rho[!is.finite(rho)]<-1
 wshape<-pmax(0,2*sqrt(matrixStats::rowVars(a)*matrixStats::rowVars(b))*(1-pmin(1,pmax(-1,rho))))
 ma<-moment_noise(a);mb<-moment_noise(b);ds<-ma$S-mb$S;dt<-ma$T-mb$T
 sr<-ds/(abs(ma$S)+abs(mb$S)+1e-8);kr<-dt/(abs(ma$T)+abs(mb$T)+1e-8)
 ws<-(.75*ds+.25*dt)/(.75*(abs(ma$S)+abs(mb$S))+.25*(abs(ma$T)+abs(mb$T))+1e-8)
 zS<-positive_ratio(ds^2,ma$vS+mb$vS);zT<-positive_ratio(dt^2,ma$vT+mb$vT)
 features<-list(KS=ne*ks^2,Wilcoxon=-log(pmax(wp,1e-300)),AD=ad,CvM=cvm,Welch=welch,Wasserstein=wass,
   Wshape=wshape,OriginalWSKR=abs(ws),SR=abs(sr),KR=abs(kr),Decancel=.75*abs(sr)+.25*abs(kr),ShapeNoise=zS+zT,SkewNoise=zS,TailNoise=zT)
 # Bounded pooled-rank basis; conditional pooled feature variances account for ties.
 u<-2*H-1;prev<-matrix(1,G,N);cur<-u;energies<-matrix(0,G,8)
 for(k in 1:8){
   if(k>1){nextp<-((2*k-1)*u*cur-(k-1)*prev)/k;prev<-cur;cur<-nextp}
   diff<-rowMeans(cur[,case,drop=FALSE])-rowMeans(cur[,control,drop=FALSE])
   energies[,k]<-positive_ratio(ne*diff^2,matrixStats::rowVars(cur))
   features[[paste0('Rank',k)]]<-energies[,k]
 }
 for(k in c(2,4,8))features[[paste0('RankEnergy',k)]]<-rowSums(energies[,1:k,drop=FALSE])
 for(q in c(.05,.1,.2,.3,.5,.7,.8,.9,.95)){
   f<-1*(H<=q);diff<-rowMeans(f[,case,drop=FALSE])-rowMeans(f[,control,drop=FALSE])
   features[[paste0('CDF',as.integer(q*100))]]<-positive_ratio(ne*diff^2,matrixStats::rowVars(f))
 }
 out<-log1p(as.matrix(as.data.frame(features)))
 out<-cbind(out,log_n_eff=rep(log(ne),G),imbalance=rep(abs(log(n1/n0)),G),tie_fraction=1-rowMeans(1/ties))
 out[!is.finite(out)]<-0;stopifnot(all(is.finite(out)));out
}
classic_features<-c('KS','Wilcoxon','AD','CvM','Welch','Wasserstein','log_n_eff','imbalance','tie_fraction')
baseline_columns<-c('OriginalWSKR','SR','KR','Decancel','ShapeNoise','KS','Wilcoxon','AD','CvM','Wasserstein','Wshape','RankEnergy2','RankEnergy4','RankEnergy8')
