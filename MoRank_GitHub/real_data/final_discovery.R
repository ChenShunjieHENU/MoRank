source('final_real_core.R')
score <- load_frozen_rule()
tcga <- tcga_data()
ref <- read_input('round7/results/all_gene_complement_audit.csv')
aux <- read_input('round13/data/TCGA_aux_Wilcoxon_rawP.csv')
locked <- discovery_ranks()[match(ref$gene,gene)]
parts <- split(seq_len(nrow(tcga$expr)),ceiling(seq_len(nrow(tcga$expr))/1500))
score_subset <- function(case,control) {
  columns <- c(control,case)
  f <- do.call(rbind,lapply(parts,function(rows) distribution_features(
    tcga$expr[rows,columns,drop=FALSE],length(control)+seq_along(case),seq_along(control))))
  score(tcga$counts[,columns,drop=FALSE],f)
}
full <- score_subset(which(tcga$group=='1'),which(tcga$group=='0'))
rank_min <- frank(-full$score,ties.method='min')
rank_max <- frank(-full$score,ties.method='max')
stopifnot(identical(as.integer(rank_max),locked$final_rank),
          cor(full$original,ref$WSKR_A_score,method='spearman')>.99999)
ref[,`:=`(MoRank_rank_min=rank_min,MoRank_rank_max=rank_max,MoRank_score=full$score,
    aux_Wilcoxon_rawP=aux$aux_Wilcoxon_rawP[match(gene,aux$gene)])]
methods <- c('MoRank','SR','DESeq2_BH','edgeR_BH')
select <- function(K) {
  possible <- sapply(methods,function(m) ref[[paste0(m,'_rank_min')]]<=K)
  definite <- sapply(methods,function(m) ref[[paste0(m,'_rank_max')]]<=K &
    switch(m,MoRank=ref$aux_Wilcoxon_rawP<.05,SR=ref$aux_Wilcoxon_rawP<.05,
           DESeq2_BH=is.finite(ref$DESeq2_q_raw),edgeR_BH=is.finite(ref$edgeR_q_raw)))
  list(definite=definite,possible=possible)
}
count_support <- with(ref,(is.finite(DESeq2_q_raw)&DESeq2_q_raw<.05) |
                       (is.finite(edgeR_q_raw)&edgeR_q_raw<.05))
budgets <- rbindlist(lapply(c(100L,250L,500L,1000L),function(K) {
  lists <- select(K)
  complement <- lists$definite[,1] & rowSums(lists$possible[,-1,drop=FALSE])==0
  data.table(K=K,MoRank_selected=sum(lists$definite[,1]),supplement=sum(complement),
             count_significant=sum(complement & count_support))
}))
lists <- select(500L)
comp <- lists$definite[,1] & rowSums(lists$possible[,-1,drop=FALSE])==0
stopifnot(identical(comp,locked$final_candidate),sum(comp)==316L)
specific <- sapply(seq_along(methods),function(i) lists$definite[,i] &
  rowSums(lists$possible[,-i,drop=FALSE])==0)
colnames(specific) <- methods
pattern <- apply(lists$definite,1,function(v) paste(as.integer(v),collapse=''))
save_result(data.table(pattern=pattern)[pattern!='0000',.N,by=pattern][order(-N,pattern)],'intersections.csv')
save_result(budgets,'budget.csv')
save_result(ref,'discovery_ranks_and_comparators.csv')
save_result(data.table(gene=ref$gene,specific,shared=lists$definite[,1]&!comp,selected=lists$definite[,1]),'membership_K500.csv')
save_result(rbindlist(lapply(budgets$K,function(K) {
  b <- select(K)
  data.table(K=K,gene=ref$gene,included=b$definite[,1]&rowSums(b$possible[,-1,drop=FALSE])==0)
})),'budget_membership.csv')
mapping <- candidate_mapping()[order(final_rank)]
x <- tcga$expr[mapping$gene,,drop=FALSE]
rownames(x) <- mapping$canonical_symbol
save_object(list(expression=x,group=tcga$group,gene_mapping=mapping),'candidate_expression.rds')
save_result(mapping,'candidate_mapping.csv')
save_result(data.table(pln_probability=full$pln_probability,mapped=nrow(mapping),
  overlap_original=sum(locked$original_candidate&comp),both_count_nonsignificant=sum(comp&
    is.finite(ref$DESeq2_q_raw)&is.finite(ref$edgeR_q_raw)&!count_support)), 'discovery_summary.csv')
cat('DISCOVERY ranks verified. Candidates',sum(comp),'mapped',nrow(mapping),'\n');flush.console()
plans <- readRDS(input('round6/raw/real/real_plans_first12.rds'))
stopifnot(length(plans)==12L)
weight <- function(s,K) {
  cut <- sort(s,decreasing=TRUE)[K]
  out <- 1*(s>cut)
  equal <- s==cut
  out[equal] <- (K-sum(out))/sum(equal)
  out
}
stability <- rbindlist(lapply(seq_along(plans),function(i) {
  plan <- plans[[i]]
  s <- score_subset(plan$case,plan$control)
  original <- readRDS(input(sprintf('round6/raw/real/real_scores_rep%d.rds',i)))
  stopifnot(identical(original$genes,ref$gene),cor(s$original,original$scores[,'WSKR_A'],method='spearman')>.99999)
  cat('STABILITY',i,'of',length(plans),'\n');flush.console()
  rbindlist(lapply(c(100L,500L,1000L),function(K) data.table(
    replicate=i,K=K,spearman=cor(s$score,full$score,method='spearman'),
    expected_overlap=sum(weight(full$score,K)*weight(s$score,K))/K,
    pln_probability=s$pln_probability,
    gene=focus_genes,rank_max=frank(-s$score,ties.method='max')[match(focus_genes,ref$gene)])))
}))
save_result(stability,'stability.csv')
