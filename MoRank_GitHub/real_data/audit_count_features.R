source("final_real_core.R")
source("../simulation/reference/count_baselines.R")
mode <- commandArgs(TRUE)[2]
if (mode=="bulk") {
 cached <- readRDS(input("round13/data/external_frozen_features.rds"))
 meta <- as.data.table(cached$metadata)
 raw <- read_input("round13/sources/GSE81089_readcounts_featurecounts.tsv.gz")
 setnames(raw,c("L608T_2122","L771T_1"),c("L608T","L771T"))
 raw <- raw[grepl("^ENSG[0-9]+$",raw[[1]])]
 sample_names <- meta$sample
 counts <- as.matrix(raw[,..sample_names]); rownames(counts) <- raw[[1]]
 storage.mode(counts) <- "integer"; rm(raw)
 map <- read_input("round13/data/ensembl_symbol_map.csv")[external_symbol_unique==TRUE]
 stopifnot(!anyDuplicated(map[!is.na(TCGA_input_gene)&nzchar(TCGA_input_gene),TCGA_input_gene]))
 indices <- match(map$ensembl[match(cached$gene,map$TCGA_input_gene)],rownames(counts))
 stopifnot(!anyNA(indices),length(indices)==21995)
 x <- log2(sweep(counts[indices,,drop=FALSE],2,colSums(counts),"/")*1e6+1)
 rownames(x) <- cached$gene
 source("../simulation/reference/redesign_core.R",chdir=TRUE)
 case <- which(meta$condition=="tumor"); control <- which(meta$condition=="normal")
 f <- distribution_features(x,case,control)
 err <- max(abs(f-cached$features))
 cat("GSE81089",nrow(counts),"denominator genes;",nrow(f),"scored; feature abs error",err,"\n"); flush.console()
 stopifnot(err<1e-10)
 old <- read_input("round13/results/external_full_cohort_ranking.csv")[match(cached$gene,gene)]
 aux <- wilcox_rows(x,case,control)
 stopifnot(max(abs(aux-old$aux_Wilcoxon_rawP))<1e-12, max(abs(f[,"SR"]-old$SR))<1e-12)
 C <- counts[indices,c(control,case),drop=FALSE]; rownames(C)<-cached$gene
 cc <- count_scores(C,length(case),length(control))
 fwrite(cbind(data.table(gene=cached$gene,aux_Wilcoxon_rawP=aux,SR=f[,"SR"]),as.data.table(cc)),
   file.path(project,"common/results/bulk_comparator_replay.csv"))
 cat("BULK_COMPARATORS errors\n");print(vapply(colnames(cc),function(n)max(abs(cc[,n]-old[[n]])),numeric(1)))
} else {
 z <- tcga_data()
 counts <- z$counts
 condition <- factor(ifelse(z$group=="1","tumor","normal"),levels=c("normal","tumor"))
 ref <- read_input("round7/results/all_gene_complement_audit.csv")
 stopifnot(identical(rownames(counts),ref$gene))
 if (mode=="tcga_deseq2") {
  dds <- DESeqDataSetFromMatrix(counts,data.frame(condition=condition,row.names=colnames(counts)),~condition)
  dds <- DESeq(dds,quiet=TRUE)
  dr <- results(dds,contrast=c("condition","tumor","normal"))
  fwrite(data.table(gene=ref$gene,q=dr$padj,logFC=dr$log2FoldChange),file.path(project,"common/results/tcga_deseq2.csv"))
 } else if (mode=="tcga_edger") {
  dge <- calcNormFactors(DGEList(counts=counts,group=condition))
  design <- model.matrix(~condition)
  dge <- estimateDisp(dge,design)
  fit <- glmQLFit(dge,design)
  tab <- topTags(glmQLFTest(fit,coef="conditiontumor"),n=Inf,sort.by="none")$table
  fwrite(data.table(gene=rownames(tab),q=tab$FDR,logFC=tab$logFC),file.path(project,"common/results/tcga_edger.csv"))
 } else if (mode=="tcga_aux") {
  source("../simulation/reference/redesign_core.R",chdir=TRUE)
  case <- which(condition=="tumor"); control <- which(condition=="normal")
  f <- distribution_features(z$expr,case,control)
  aux <- exp(-expm1(f[,"Wilcoxon"]))
  raw <- read_input("round13/data/TCGA_aux_Wilcoxon_rawP.csv")
  stopifnot(identical(raw$gene,ref$gene),max(abs(raw$aux_Wilcoxon_rawP-aux))<1e-12)
  sr <- read_input("round3/inputs/gene_scores.csv")[match(ref$gene,gene)]
  stopifnot(max(abs(f[,"SR"]-sr$abs_SR))<1e-10)
  cat("TCGA_AUX_SR exact ranking inputs",length(aux),"\n")
 } else stop(mode)
}

