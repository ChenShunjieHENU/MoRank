source('final_real_core.R')
suppressPackageStartupMessages(library(edgeR))
source('functional_core.R')
candidate <- candidate_mapping()[order(final_rank)]
rawmap <- read_input('round10/data/raw_gene_mapping.csv')[canonical_collision==FALSE]
stopifnot(!anyDuplicated(rawmap$canonical_symbol))
candidate[,raw_symbol:=rawmap$raw_symbol[match(canonical_symbol,rawmap$canonical_symbol)]]
display <- head(unique(c(focus_genes,candidate[!is.na(raw_symbol),canonical_symbol])),14L)
epi <- readRDS(input('round10/data/classified_epithelial_counts.rds'))
pb <- readRDS(input('round10/data/pseudobulk_full.rds'))
tab <- as.data.table(pb$metadata)
pairs <- sort(unique(tab$patient))
stopifnot(length(pairs)==9L,all(table(tab$patient)==2L))
normal <- match(paste(pairs,'normal'),paste(tab$patient,tab$condition))
tumor <- match(paste(pairs,'tumor'),paste(tab$patient,tab$condition))
delta <- function(m) m[,tumor,drop=FALSE]-m[,normal,drop=FALSE]
genes <- candidate[!is.na(raw_symbol)&raw_symbol %in% rownames(epi$counts)]
group_matrix <- function(meta,sample_ids) sparseMatrix(i=seq_len(nrow(meta)),
    j=match(meta$Sample,sample_ids),x=1,dims=c(nrow(meta),length(sample_ids)))
cols <- which(epi$metadata$Sample %in% tab$Sample)
meta <- epi$metadata[cols]
Y <- epi$counts[genes$raw_symbol,cols,drop=FALSE]
rownames(Y) <- genes$gene
G <- group_matrix(meta,tab$Sample)
stopifnot(all(colSums(G)>=20),all(as.numeric(colSums(epi$counts[,cols,drop=FALSE]))==meta$library_UMI))
Z <- Y %*% Diagonal(x=10000/meta$library_UMI)
Z@x <- log1p(Z@x)
observed <- sweep(as.matrix((Y>0)%*%G),2,colSums(G),'/')
mean_exp <- sweep(as.matrix(Z%*%G),2,colSums(G),'/')
detection <- data.table(gene=genes$gene,endpoint='detection',paired_stats(delta(observed)))
expression <- data.table(gene=genes$gene,endpoint='mean_logCP10K',paired_stats(delta(mean_exp)))

fixed <- which(meta$library_UMI>=1000L)
GD <- group_matrix(meta[fixed],tab$Sample)
stopifnot(all(colSums(GD)>=20))
raw_detection <- as.matrix(Y[,fixed,drop=FALSE])
L <- meta$library_UMI[fixed]
expected <- t(vapply(seq_len(nrow(raw_detection)),function(i) {
  y <- raw_detection[i,]
  phyper(0,y,L-y,1000L,lower.tail=FALSE)
},numeric(length(fixed))))
stopifnot(all(is.finite(expected)),all(expected>=0 & expected<=1))
depth <- sweep(as.matrix(expected%*%GD),2,colSums(GD),'/')
ds <- data.table(gene=genes$gene,paired_stats(delta(depth)))
ds <- merge(candidate_list()[,.(gene)],ds,by='gene',all.x=TRUE)
ds[,q_candidate:=finite_bh(p)]
save_result(ds,'epithelial_depth1000.csv')
row_irs <- match('IRS2',genes$gene)
save_result(rbind(data.table(tab[,.(patient,condition)],detection=observed[row_irs,],scale='Observed'),
  data.table(tab[,.(patient,condition)],detection=depth[row_irs,],scale='Fixed 1,000 UMI')),'IRS2_depth_pairs.csv')

design <- model.matrix(~factor(patient)+factor(condition,levels=c('normal','tumor')),tab)
dge <- DGEList(pb$counts)
dge <- dge[filterByExpr(dge,design=design),,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge,design,robust=TRUE)
fit <- glmQLFit(dge,design,robust=TRUE)
e <- as.data.table(topTags(glmQLFTest(fit,coef=ncol(design)),n=Inf,sort.by='none')$table,keep.rownames='raw_symbol')
audit <- read_input('round10/results/pseudobulk_all_tested_genes.csv')
stopifnot(nrow(e)==nrow(audit),max(abs(e$PValue-audit$PValue[match(e$raw_symbol,audit$raw_symbol)]))<1e-10)
selected <- merge(candidate_list()[,.(gene)],candidate[,.(gene,raw_symbol)],by='gene',all.x=TRUE)
selected <- merge(selected,e,by='raw_symbol',all.x=TRUE)
endpoints <- merge(CJ(gene=candidate_list()$gene,endpoint=c('detection','mean_logCP10K','pseudobulk_edgeR')),
  rbind(detection[,.(gene,endpoint,p)],expression[,.(gene,endpoint,p)],
        selected[,.(gene,endpoint='pseudobulk_edgeR',p=PValue)]),by=c('gene','endpoint'),all.x=TRUE)
endpoints[,q_joint:=finite_bh(p)]
stopifnot(nrow(endpoints)==948L)
save_result(endpoints,'epithelial_three_endpoints.csv')
selected[,q_joint:=endpoints[endpoint=='pseudobulk_edgeR',q_joint[match(selected$gene,gene)]]]
save_result(selected,'epithelial_pseudobulk.csv')
save_result(merge(rbind(detection,expression),endpoints[,.(gene,endpoint,q_joint)],by=c('gene','endpoint')),
  'epithelial_paired_effects.csv')
irs <- endpoints[gene=='IRS2']
irs[,holm:=p.adjust(ifelse(is.finite(p),p,1),'holm')]
save_result(irs,'IRS2_endpoint_tests.csv')

# Recompute patient-weighted lineage expression for the revised candidate display.
rm(epi,Y,Z,expected,raw_detection)
invisible(gc())
cache <- readRDS(input('round10/data/primary_lung_input_cache.rds'))
cell <- read_input('round10/data/primary_lung_cell_metadata.csv.gz')
stopifnot(identical(colnames(cache$counts),cell$Index),nrow(cell)==88144L)
sub <- which(cell$condition=='tumor' & !cell$cell_group %in% c('Undetermined','Unassigned epithelial'))
cell <- cell[sub]
raw_symbols <- rawmap$raw_symbol[match(display,rawmap$canonical_symbol)]
present <- !is.na(raw_symbols)&raw_symbols %in% rownames(cache$counts)
counts <- cache$counts[raw_symbols[present],sub,drop=FALSE]
rownames(counts) <- display[present]
rm(cache)
norm <- counts%*%Diagonal(x=10000/cell$library_UMI)
norm@x <- log1p(norm@x)
groups <- unique(cell[,.(patient,cell_group)])
id <- match(paste(cell$patient,cell$cell_group),paste(groups$patient,groups$cell_group))
group <- sparseMatrix(i=seq_len(nrow(cell)),j=id,x=1,dims=c(nrow(cell),nrow(groups)))
sizes <- colSums(group)
D <- sweep(as.matrix((counts>0)%*%group),2,sizes,'/')
E <- sweep(as.matrix(norm%*%group),2,sizes,'/')
loc <- rbindlist(lapply(which(sizes>=20),function(i) data.table(gene=rownames(D),
  patient=groups$patient[i],cell_group=groups$cell_group[i],det=D[,i],expr=E[,i])))
loc <- loc[,.(patients=uniqueN(patient),patient_mean_detection=mean(det),
  patient_mean_logCP10K=mean(expr)),by=.(gene,cell_group)][patients>=3]
loc[,expression_z:=if(sd(patient_mean_logCP10K)>0) as.numeric(scale(patient_mean_logCP10K)) else 0,by=gene]
save_result(loc,'epithelial_lineage_display.csv')
save_result(data.table(mapped=nrow(genes),pseudobulk_tested=sum(is.finite(selected$PValue)),
  pseudobulk_joint_BH=sum(selected$q_joint<.05),depth_tested=sum(is.finite(ds$p)),
  depth_BH=sum(ds$q_candidate<.05),patients=length(pairs),fixed_depth_cells=length(fixed)),'epithelial_summary.csv')
cat('Single-cell paired and lineage results complete\n')
