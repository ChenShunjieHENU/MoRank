# Whole-pathway patient activity, not candidate-only scores and not causal perturbations.
pathrows<-rbind(tg[database=='Hallmark'],tg[database=='KEGG' & term_id %in% kctx]);paths<-split(pathrows$gene,pathrows$term_id);paths<-lapply(paths,unique);stopifnot(length(paths)==75)
activity<-list();patientvals<-list()
score_cohort<-function(counts,meta,symbols,cohort){
 stopifnot(length(symbols)==nrow(counts));ok<-!is.na(symbols)&nzchar(symbols)&!duplicated(symbols)&!duplicated(symbols,fromLast=TRUE);counts<-counts[ok,,drop=FALSE];rownames(counts)<-symbols[ok]
 # Full measured unambiguous genes form the rank denominator in each sample.
 xx<-log2(sweep(counts,2,colSums(counts),'/')*1e6+1);Q<-apply(xx,2,rank,ties.method='average')/(nrow(xx)+1);Q<-Q-rowMeans(Q)
 ni<-which(meta$condition=='normal');ti<-match(as.character(meta$patient[ni]),as.character(meta$patient[meta$condition=='tumor']));ti<-which(meta$condition=='tumor')[ti];stopifnot(length(ni)==9,!anyNA(ti))
 S<-t(vapply(paths,function(gs){ix<-which(rownames(Q) %in% gs);if(length(ix)<5)rep(NA_real_,ncol(Q)) else colMeans(Q[ix,,drop=FALSE])},numeric(ncol(Q))));rownames(S)<-names(paths);colnames(S)<-meta$Sample
 delta<-S[,ti,drop=FALSE]-S[,ni,drop=FALSE];st<-as.data.table(paired_stats(delta,B=5000L,seed=202611083));st[,`:=`(term_id=names(paths),cohort=cohort,measured_genes=vapply(paths,function(gs)sum(rownames(Q) %in% gs),integer(1)))];st[,q:=p.adjust(ifelse(is.finite(p),p,1),'BH')]
 val<-data.table(term_id=rep(names(paths),each=ncol(S)),sample=rep(meta$Sample,times=nrow(S)),score=as.vector(t(S)));val<-merge(val,meta[,.(sample=Sample,patient=as.character(patient),condition)],by='sample');val[,cohort:=cohort]
 list(stats=st,values=val,log_expression=xx,metadata=meta)
}
sc<-readRDS(file.path(P,'round10/data/pseudobulk_full.rds'));sg<-rd('round10/data/raw_gene_mapping.csv');smatch<-sg$canonical_symbol[match(rownames(sc$counts),sg$raw_symbol)];sc$metadata[,condition:=as.character(condition)];a<-score_cohort(sc$counts,sc$metadata,smatch,'Epithelium');activity[[1]]<-a$stats;patientvals[[1]]<-a$values
# Preserve actual paired epithelial display values for all five fixed examples.
vv<-as.data.table(as.table(a$log_expression[intersect(c('ABCC3','NQO1','IRS2','KIT','FGFR3'),rownames(a$log_expression)),,drop=FALSE]));setnames(vv,c('gene','sample','logCPM_plus1'));vv<-merge(vv,a$metadata[,.(sample=Sample,patient=as.character(patient),condition)],by='sample');out(vv,'epithelial_paired_gene_values.csv')
meta<-as.data.table(readRDS(file.path(P,'round13/data/external_analysis_manifest.rds'))$paired_metadata);setnames(meta,'sample','Sample');raw<-rd('round13/sources/GSE81089_readcounts_featurecounts.tsv.gz');setnames(raw,c('L608T_2122','L771T_1'),c('L608T','L771T'));raw<-raw[grepl('^ENSG[0-9]+$',raw[[1]])];sns<-meta$Sample;counts<-as.matrix(raw[,..sns]);rownames(counts)<-raw[[1]];mm<-rd('round13/data/ensembl_symbol_map.csv');symbols<-mm$canonical_symbol[match(rownames(counts),mm$ensembl)];b<-score_cohort(counts,meta,symbols,'Tissue');activity[[2]]<-b$stats;patientvals[[2]]<-b$values
act<-rbindlist(activity);nm<-unique(pathrows[,.(term_id,term_name)]);act<-merge(act,nm,by='term_id',all.x=TRUE);out(act,'pathway_activity_tests.csv');out(rbindlist(patientvals),'pathway_patient_scores.csv')
