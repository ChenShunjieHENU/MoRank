suppressPackageStartupMessages({library(DESeq2); library(edgeR)})
count_scores<-function(counts,n1,n0){
 group<-factor(c(rep('control',n0),rep('case',n1)),levels=c('control','case'))
 meta<-data.frame(condition=group,row.names=colnames(counts))
 dds<-DESeqDataSetFromMatrix(round(counts),meta,~condition)
 dds<-DESeq(dds,quiet=TRUE,parallel=FALSE)
 dr<-as.data.frame(results(dds,contrast=c('condition','case','control')))
 dge<-DGEList(counts=round(counts),group=group);keep<-filterByExpr(dge,group=group)
 stopifnot(sum(keep)>10);dge<-dge[keep,,keep.lib.sizes=FALSE];dge<-calcNormFactors(dge);design<-model.matrix(~group)
 dge<-estimateDisp(dge,design);fit<-glmQLFit(dge,design);et<-topTags(glmQLFTest(fit,coef=2),n=Inf,sort.by='none')$table
 sc<-function(p){p[!is.finite(p)]<-1;-log10(pmax(p,1e-300))}
 ep<-eq<-rep(1,nrow(counts));ep[keep]<-et$PValue;eq[keep]<-et$FDR
 cbind(DESeq2_P=sc(dr$pvalue),DESeq2_BH=sc(dr$padj),edgeR_P=sc(ep),edgeR_BH=sc(eq))
}
