suppressPackageStartupMessages(library(matrixStats));source('uniform_permutation.R');source('redesign_core.R')
set.seed(117);state<-.Random.seed
u<-rbind(c(.1,.1,.3,.4),c(.4,.3,.2,.1));z<-priority_permutations(u,1L,1L)
stopifnot(identical(z$tied_rows,1L),all(apply(z$indices,1,function(v)identical(sort(as.integer(v)),1:4))),identical(state,.Random.seed))
a<-distribution_features(matrix(c(0,0,0,1),1),3:4,1:2)
stopifnot(abs(expm1(a[1,'AD'])-1)<1e-12,abs(expm1(a[1,'CvM'])-1/16)<1e-12,abs(expm1(a[1,'KS'])-.25)<1e-12)
cat('PASS forced-priority-tie native-uniform fallback,permutation membership,RNGstream preservation,hand-calculated discreteAD/CvM/KS\n')
