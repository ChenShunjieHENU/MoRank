# Distinct random priorities give a uniform permutation. A finite-RNG priority tie is replaced
# by an independently seeded native sample.int permutation; the priority RNG stream is preserved.
priority_permutations<-function(u,cell,b){
 N<-ncol(u);r<-matrixStats::rowRanks(u,ties.method='max',useNames=FALSE)
 bad<-which(rowSums(r)!=N*(N+1)/2)
 if(length(bad)){
   saved_seed<-.Random.seed
   for(i in bad){set.seed(900000000L+cell*1000000L+b*1000L+i);r[i,]<-sample.int(N,N,replace=FALSE)}
   .Random.seed<<-saved_seed
 }
 list(indices=r,tied_rows=bad)
}
