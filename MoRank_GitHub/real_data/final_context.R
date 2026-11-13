source('final_real_core.R')
suppressPackageStartupMessages({library(igraph); library(jsonlite); library(edgeR); library(matrixStats)})
source('functional_core.R')
set.seed(202611084L)
cand <- candidate_mapping()$canonical_symbol
sig <- read_result('bulk_paired_candidates.csv')[q_candidate<.05 & is.finite(PValue),canonical_symbol]
sig <- unique(na.omit(sig))
mapped <- read_input('round9/results/mapped_universe.csv')
terms <- read_input('round9/results/term_gene.csv')[database %in% c('GO_BP','Reactome')]
backgrounds <- read_input('round9/results/database_backgrounds.csv')
tg <- unique(read_input('round20_network/data/full_term_gene.csv'))
plan <- read_json('pathway_plan.json')
kctx <- unlist(plan$context_kegg)
hctx <- unlist(plan$context_hallmark)

bin <- function(v) {
  cuts <- unique(quantile(v,seq(0,1,.2),names=FALSE))
  if(length(cuts)==1L) return(rep(1L,length(v)))
  findInterval(v,cuts,all.inside=TRUE,rightmost.closed=TRUE)
}
ora <- list()
matched <- list()
for(db in c('GO_BP','Reactome','KEGG','Hallmark')) {
  context <- if(db=='KEGG') kctx else if(db=='Hallmark') hctx else character()
  if(db %in% c('GO_BP','Reactome')) {
    t <- unique(terms[database==db])
    u <- sort(unique(backgrounds[database==db,gene]))
  } else {
    t <- tg[database==db & gene %in% mapped$gene]
    u <- sort(unique(t$gene))
    eligible <- t[,.(n=uniqueN(gene)),by=term_id][n>=10 & (db=='Hallmark' | n<=500),term_id]
    t <- t[term_id %in% eligible]
  }
  names <- unique(t[,.(term_id,term_name)])[order(term_id)]
  A <- sparseMatrix(i=match(t$term_id,names$term_id),j=match(t$gene,u),x=1,dims=c(nrow(names),length(u)))
  K <- as.integer(rowSums(A))
  one <- function(genes,label) {
    use <- which(u %in% genes)
    n <- length(use)
    hits <- as.integer(rowSums(A[,use,drop=FALSE]))
    annotation <- t[gene %in% genes,.(hit_genes=paste(sort(unique(gene)),collapse=';')),by=term_id]
    m <- data.table(names,database=db,set=label,background_n=length(u),K=K,set_annotated_n=n,
        hit_n=hits,p=phyper(hits-1,K,length(u)-K,n,lower.tail=FALSE),
        fold=hits/n/(K/length(u)),context=names$term_id %in% context)
    m <- merge(m,annotation,by='term_id',all.x=TRUE,sort=FALSE)
    m[is.na(hit_genes),hit_genes:='']
    m
  }
  rr <- rbind(one(cand,'Final candidates'),one(sig,'Bulk supported'))
  rr[,q:=p.adjust(p,'BH')]
  ora[[db]] <- rr
  cv <- mapped[match(u,gene)]
  stopifnot(identical(cv$gene,u))
  strata <- paste(bin(cv$mean_logCPM),bin(cv$detection_logCPM_gt0),sep=':')
  selected <- which(u %in% cand)
  nt <- table(strata[selected]); pools <- split(seq_along(u),strata)
  B <- if(db=='Hallmark') 10000L else 1000L
  set.seed(202611084L)
  draws <- replicate(B,unlist(lapply(names(nt),function(key) {
    pool <- pools[[key]]
    pool[sample.int(length(pool),nt[[key]])]
  }),use.names=FALSE),simplify=FALSE)
  D <- sparseMatrix(i=unlist(draws,use.names=FALSE),j=rep(seq_len(B),lengths(draws)),x=1,dims=c(length(u),B))
  H <- as.matrix(A%*%D)
  obs <- as.integer(rowSums(A[,selected,drop=FALSE]))
  quant <- t(apply(H,1,quantile,c(.025,.5,.975),names=FALSE))
  tail <- rowSums(H>=obs)
  m <- data.table(names,database=db,hit_n=obs,reference_mean=rowMeans(H),lo=quant[,1],median=quant[,2],
                  hi=quant[,3],exceedances=tail,p=(1+tail)/(B+1),draws=B)
  m[,q:=p.adjust(p,'BH')]
  matched[[db]] <- m
  if(db=='Hallmark') {
    display <- c('HALLMARK_XENOBIOTIC_METABOLISM','HALLMARK_TNFA_SIGNALING_VIA_NFKB')
    save_result(rbindlist(lapply(display,function(term) data.table(term_id=term,
      hits=H[match(term,names$term_id),]))),'hallmark_matched_draws.csv')
  }
  cat('Enrichment',db,nrow(rr),'tests',sum(rr$q<.05),'BH results\n');flush.console()
}
save_result(rbindlist(ora),'enrichment_all.csv')
save_result(rbindlist(matched),'enrichment_matched.csv')

# Re-run the candidate-induced correlation graphs and the fixed sample-size sensitivity.
tcga <- load_object('candidate_expression.rds')
X <- tcga$expression
pairs <- which(upper.tri(matrix(0,nrow(X),nrow(X))),arr.ind=TRUE)
edge_stats <- function(x) {
  n <- ncol(x)
  rho <- suppressWarnings(cor(t(x),method='spearman'))[pairs]
  p <- 2*pt(abs(rho)*sqrt((n-2)/pmax(1-rho^2,1e-15)),df=n-2,lower.tail=FALSE)
  q <- finite_bh(p)
  data.table(from=rownames(x)[pairs[,1]],to=rownames(x)[pairs[,2]],rho=rho,q=q,
             keep=is.finite(rho)&abs(rho)>=.5&q<.01)
}
x_tumor <- X[,tcga$group=='1',drop=FALSE]
tumor <- edge_stats(x_tumor)
normal <- edge_stats(X[,tcga$group=='0',drop=FALSE])
vertices <- data.table(name=rownames(X))
graphs <- list(Tumor=graph_from_data_frame(tumor[keep==TRUE,.(from,to,weight=abs(rho))],FALSE,vertices),
               Normal=graph_from_data_frame(normal[keep==TRUE,.(from,to,weight=abs(rho))],FALSE,vertices))
save_object(graphs,'coexpression_graphs.rds')
save_result(merge(tumor[,.(from,to,rho_tumor=rho,keep_tumor=keep)],
    normal[,.(from,to,rho_normal=rho,keep_normal=keep)],by=c('from','to')),'coexpression_pairs.csv')
set.seed(202611084L)
save_result(rbindlist(lapply(seq_len(100L),function(i) data.table(replicate=i,
   edges=sum(edge_stats(x_tumor[,sample.int(ncol(x_tumor),59L),drop=FALSE])$keep)))), 'coexpression_resamples.csv')

mapped_string <- read_result('string_mapped.csv')
networks <- lapply(c('functional','physical'),function(kind) {
  ed <- read_result(paste0('string_',kind,'.csv'))[score>=.7,.(from,to,weight=score)]
  graph_from_data_frame(ed,FALSE,data.table(name=mapped_string$gene))
})
names(networks) <- c('functional','physical')
save_object(networks,'string_graphs.rds')
save_result(rbindlist(lapply(names(networks),function(kind) {
  g <- networks[[kind]]
  data.table(gene=V(g)$name,type=kind,degree=as.integer(degree(g)),
             strength=as.numeric(strength(g)),external_support=V(g)$name %in% sig)
})),'string_node_metrics.csv')
save_result(data.table(mapped=nrow(X),pairs=nrow(tumor),string_mapped=nrow(mapped_string),
  functional=ecount(networks$functional),physical=ecount(networks$physical),
  tumor=ecount(graphs$Tumor),normal=ecount(graphs$Normal)),'network_summary.csv')

# These scores use all measured members, so the tested pathways and input cohorts stay fixed.
P <- parent
rd <- read_input
out <- save_result
source('pathway_activity.R')

# The same paired participants and IRS2 are retained, so rerun the six fixed module tests.
module_pairs <- read_input('round10/results/module_IRS2_paired_changes.csv')
module_tests <- module_pairs[,{
  shift <- paired_stats(matrix(delta_module,nrow=1))
  association <- spearman_perm(IRS2_delta_mean_logCP10K,delta_module)
  .(test=c('paired_module_shift','IRS2_change_correlation'),
    estimate=c(shift$mean_delta,association[['rho']]),p=c(shift$p,association[['p']]),
    ci_low=c(shift$ci_low,NA_real_),ci_high=c(shift$ci_high,NA_real_),n_pairs=.N)
},by=module]
module_tests[,q_BH6:=finite_bh(p)]
reference_modules <- read_input('round10/results/module_tests_all6.csv')
check <- merge(module_tests,reference_modules,by=c('module','test'))
stopifnot(nrow(check)==6L,max(abs(check$q_BH6.x-check$q_BH6.y))<1e-12)
save_result(module_tests,'IRS2_module_tests.csv')

# IRS2 is still selected. Recalculate the original matched perturbation endpoints
# directly from the frozen target/sham/control response vectors and the fixed modules.
stopifnot('IRS2' %in% cand)
ko <- readRDS(input('round10/results/reference_KO/perturbations_all102.rds'))
sub <- readRDS(input('round10/data/IRS2_reference_subnetwork.rds'))
stopifnot(length(ko)==102L,length(sub$controls)==100L,length(sub$response_genes)==1899L)
endpoints <- c(list(global=sub$response_genes),sub$modules)
values <- rbindlist(lapply(names(ko),function(target) {
  d <- ko[[target]]$distances
  v <- vapply(endpoints,function(genes) mean(d[genes]),numeric(1))
  v['global'] <- sum(d[sub$response_genes])
  stopifnot(max(abs(v-ko[[target]]$statistics[names(v)]))<1e-12)
  data.table(target=target,endpoint=names(v),value=as.numeric(v))
}))
summary <- values[,{
  v <- value[target %in% sub$controls]
  obs <- value[target=='IRS2']
  .(observed=obs,sham=value[target=='SHAM'],median=median(v),lo=quantile(v,.025),
    hi=quantile(v,.975),p=(1+sum(v>=obs))/(1+length(v)))
},by=endpoint]
summary[,holm:=p.adjust(p,'holm')]
save_result(values,'IRS2_reference_values.csv')
save_result(summary,'IRS2_reference_summary.csv')
print(summary)
