suppressPackageStartupMessages(library(matrixStats))
moments <- function(x,eps=1e-8) {
  mu<-rowMeans(x);z<-sweep(x,1,mu,'-');sdv<-sqrt(rowMeans(z*z))
  bad<-!is.finite(sdv)|sdv<=eps;sdv[bad]<-1
  z<-sweep(z,1,sdv,'/');z2<-z*z
  skew<-rowMeans(z2*z);kurt<-rowMeans(z2*z2)-3
  tail<-sign(kurt)*log1p(abs(kurt));skew[bad]<-0;tail[bad]<-0
  list(skew=skew,tail=tail,constant=bad)
}
shape_scores <- function(case,control,alpha=.75,eps=1e-8) {
  a<-moments(case,eps);b<-moments(control,eps)
  ds<-a$skew-b$skew;dt<-a$tail-b$tail
  A<-abs(a$skew)+abs(b$skew);B<-abs(a$tail)+abs(b$tail)
  sr<-ds/(A+eps);kr<-dt/(B+eps)
  list(WSKR=(alpha*ds+(1-alpha)*dt)/(alpha*A+(1-alpha)*B+eps),SR=sr,KR=kr,
       effective_skew_weight=alpha*(A+eps)/(alpha*A+(1-alpha)*B+eps),
       signed_simple=alpha*sr+(1-alpha)*kr,
       magnitude_simple=alpha*abs(sr)+(1-alpha)*abs(kr))
}
wilcox_rows <- function(x,case_cols,control_cols) {
  stopifnot(length(case_cols)+length(control_cols)==ncol(x),length(intersect(case_cols,control_cols))==0)
  n1<-length(case_cols);n2<-length(control_cols);N<-n1+n2
  ranks<-matrixStats::rowRanks(x,ties.method='average',useNames=FALSE)
  U<-rowSums(ranks[,case_cols,drop=FALSE])-n1*(n1+1)/2
  # The deficit in squared ranks gives the exact pooled tie correction.
  tie_sum<-12*(N*(N+1)*(2*N+1)/6-rowSums(ranks*ranks))
  variance<-n1*n2/12*((N+1)-tie_sum/(N*(N-1)))
  centered<-U-n1*n2/2
  z<-(centered-sign(centered)*.5)/sqrt(pmax(variance,0))
  p<-2*pnorm(-abs(z));p[!is.finite(p)|variance<=0]<-1
  p
}
auc <- function(score,truth) {
  stopifnot(length(score)==length(truth),all(is.finite(score)))
  m<-sum(truth);n<-sum(!truth)
  if(m==0||n==0)return(NA_real_)
  (sum(rank(score,ties.method='average')[truth])-m*(m+1)/2)/(m*n)
}
