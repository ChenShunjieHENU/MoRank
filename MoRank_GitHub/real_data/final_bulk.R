source('final_real_core.R')
suppressPackageStartupMessages(library(edgeR))
score <- load_frozen_rule()
cached <- readRDS(input('round13/data/external_frozen_features.rds'))
meta <- as.data.table(cached$metadata)
raw <- read_input('round13/sources/GSE81089_readcounts_featurecounts.tsv.gz')
setnames(raw,c('L608T_2122','L771T_1'),c('L608T','L771T'))
raw <- raw[grepl('^ENSG[0-9]+$',raw[[1]])]
sample_names <- meta$sample
counts <- as.matrix(raw[,..sample_names])
rownames(counts) <- raw[[1]]
stopifnot(ncol(counts)==127L,length(cached$gene)==21995L)
rm(raw)
old <- read_input('round13/results/external_full_cohort_ranking.csv')[match(cached$gene,gene)]
prediction <- score(counts,cached$features)
stopifnot(cor(prediction$original,old$WSKR_A,method='spearman')>.99999)
old[,`:=`(MoRank_score=prediction$score,MoRank_rank_min=frank(-prediction$score,ties.method='min'),
          MoRank_rank_max=frank(-prediction$score,ties.method='max'))]
old[,complement:=MoRank_rank_max<=500 & aux_Wilcoxon_rawP<.05 &
      pmin(SR_rank_min,DESeq2_BH_rank_min,edgeR_BH_rank_min)>500]
save_result(old[,.(gene,MoRank_score,MoRank_rank_min,MoRank_rank_max,complement)],'GSE81089_reranking.csv')
cand <- candidate_list()
save_result(data.table(external_candidates=sum(old$complement),overlap=sum(old[complement==TRUE,gene] %in% cand$gene),
    pln_probability=prediction$pln_probability), 'external_reranking_summary.csv')

pairids <- intersect(meta[condition=='normal',patient],meta[condition=='tumor',patient])
pm <- meta[patient %in% pairids][order(patient,condition)]
pm[,`:=`(patient=factor(patient),condition=factor(condition,levels=c('normal','tumor')))]
stopifnot(length(pairids)==9L,all(table(pm$patient)==2L))
design <- model.matrix(~patient+condition,pm)
dge <- DGEList(counts[,pm$sample,drop=FALSE])
dge <- dge[filterByExpr(dge,design=design),,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge,design,robust=TRUE)
fit <- glmQLFit(dge,design,robust=TRUE)
tab <- as.data.table(topTags(glmQLFTest(fit,coef='conditiontumor'),n=Inf,sort.by='none')$table,keep.rownames='ensembl')
mapping <- read_input('round13/data/ensembl_symbol_map.csv')[nzchar(TCGA_input_gene)&!is.na(TCGA_input_gene)&
    external_symbol_unique==TRUE & ensembl %in% rownames(counts)]
stopifnot(!anyDuplicated(mapping$TCGA_input_gene),nrow(tab)==24883L)
audit <- read_input('round13/results/external_paired_all_tested.csv.gz')
stopifnot(max(abs(tab$PValue-audit$PValue[match(tab$ensembl,audit$ensembl)]))<1e-10)
attach_tests <- function(genes) {
  m <- mapping[match(genes,TCGA_input_gene),.(ensembl,canonical_symbol)]
  t <- tab[match(m$ensembl,ensembl),.(logFC,logCPM,F,PValue,FDR)]
  cbind(data.table(gene=genes),m,t)
}
ref <- read_input('round7/results/all_gene_complement_audit.csv')
locked <- discovery_ranks()
out <- attach_tests(cand$gene)
out[,`:=`(q_candidate=finite_bh(PValue),discovery_logFC=ref$DESeq2_log2FC[match(gene,ref$gene)])]
out[,both_count_nonsignificant:=ref[match(out$gene,ref$gene),
    is.finite(DESeq2_q_raw)&DESeq2_q_raw>=.05&is.finite(edgeR_q_raw)&edgeR_q_raw>=.05]]
out[,direction:=fcase(!is.finite(PValue),'Unavailable',q_candidate>=.05,'NS',
    sign(logFC)==sign(discovery_logFC),'Concordant',default='Discordant')]
save_result(out,'bulk_paired_candidates.csv')
groups <- list('Final MoRank, cap 500'=cand$gene,
  'Final MoRank, matched budget'=head(cand$gene,261L),
  'Original MoRank, cap 500'=locked[original_candidate==TRUE,gene])
summary <- rbindlist(lapply(names(groups),function(name) {
  m <- attach_tests(groups[[name]])
  m[,q:=finite_bh(PValue)]
  strict <- ref[match(m$gene,ref$gene),is.finite(DESeq2_q_raw)&DESeq2_q_raw>=.05&
    is.finite(edgeR_q_raw)&edgeR_q_raw>=.05]
  m[,.(set=name,n=.N,mapped=sum(!is.na(ensembl)),tested=sum(is.finite(PValue)),
       bh=sum(q<.05),fdr=sum(FDR<.05,na.rm=TRUE),higher=sum(q<.05&logFC>0,na.rm=TRUE),
       lower=sum(q<.05&logFC<0,na.rm=TRUE),strict=sum(strict),strict_supported=sum(strict&q<.05))]
}))
stopifnot(identical(summary$bh,c(175L,139L,79L)))
save_result(summary,'bulk_paired_summary.csv')
lcpm <- cpm(dge,log=TRUE,prior.count=2)
vals <- out[is.finite(PValue),.(gene,ensembl)]
save_result(rbindlist(lapply(seq_len(nrow(vals)),function(j) data.table(
  gene=vals$gene[j],patient=as.character(pm$patient),condition=as.character(pm$condition),
  logCPM=as.numeric(lcpm[vals$ensembl[j],])))),'bulk_paired_values.csv')
print(summary)
