suppressPackageStartupMessages({library(data.table); library(Matrix)})
args <- commandArgs(TRUE)
root <- args[1]
src <- file.path(Sys.getenv("MORANK_INPUT_DIR"), "round10")
source_counts <- file.path(src,"sources/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
X <- readRDS(source_counts)
if (is.data.frame(X)) {
 bad <- which(!vapply(X,is.numeric,logical(1)))
 if (identical(bad, 1L)) { ids <- X[[1]]; X[[1]] <- NULL; rownames(X) <- ids }
 else stopifnot(!length(bad))
 X <- as.matrix(X)
}
if (!inherits(X,"dgCMatrix")) X <- as(X,"dgCMatrix")
meta <- fread(file.path(src,"sources/cell_annotation.tsv"), na.strings=NULL)
stopifnot(!anyDuplicated(meta$Index),setequal(meta$Index,colnames(X)))
meta <- meta[match(colnames(X), Index)]
rawdim <- dim(X)
i <- which(meta$Sample_Origin %in% c("nLung","tLung"))
X <- X[,i,drop=FALSE]; meta <- meta[i]
saved <- readRDS(file.path(src,"data/primary_lung_input_cache.rds"))
equal_sparse <- function(a,b) identical(dimnames(a),dimnames(b)) && identical(a@i,b@i) &&
 identical(a@p,b@p) && identical(a@x,b@x)
stopifnot(identical(saved$original_dimension,rawdim), equal_sparse(saved$counts,X),
 identical(as.data.frame(meta),as.data.frame(saved$metadata)))
rm(saved); gc()
meta[,patient:=sprintf("P%04d",as.integer(sub("^LUNG_[NT]","",Sample)))]
meta[,condition:=ifelse(Sample_Origin=="tLung","tumor","normal")]
meta[,valid_epi:=Cell_type=="Epithelial cells" &
 ((condition=="normal" & Cell_subtype %in% c("AT1","AT2","Ciliated","Club")) |
  (condition=="tumor" & Cell_subtype %in% c("tS1","tS2","tS3","Malignant cells")))]
epi <- readRDS(file.path(src,"data/classified_epithelial_counts.rds"))
stopifnot(equal_sparse(X[,meta$valid_epi,drop=FALSE],epi$counts),
 identical(meta[valid_epi==TRUE,Index],epi$metadata$Index))
by_patient <- dcast(meta[valid_epi==TRUE,.(n_cells=.N),by=.(patient,condition)],
 patient~condition,value.var="n_cells",fill=0)
pairs <- sort(by_patient[normal>=20 & tumor>=20,patient])
eligible <- which(meta$valid_epi & meta$patient %in% pairs)
selected <- meta[eligible]
sampletab <- unique(selected[,.(patient,Sample,condition)])
setorder(sampletab,patient,condition)
G <- sparseMatrix(i=seq_len(nrow(selected)),j=match(selected$Sample,sampletab$Sample),
 x=1,dims=c(nrow(selected),nrow(sampletab)))
sampletab[,n_cells:=as.numeric(colSums(G))]
PB <- as.matrix(X[,eligible,drop=FALSE] %*% G)
colnames(PB) <- sampletab$Sample
saved_pb <- readRDS(file.path(src,"data/pseudobulk_full.rds"))
stopifnot(identical(PB,saved_pb$counts),identical(as.data.frame(sampletab),as.data.frame(saved_pb$metadata)))
cat("RAW_SPARSE_COUNTS_EXACT",rawdim,"primary lung",ncol(X),"classified epithelial",
 sum(meta$valid_epi),"paired patients",length(pairs),"pseudobulk",dim(PB),"\n")
