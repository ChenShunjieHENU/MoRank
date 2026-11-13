project <- normalizePath(commandArgs(TRUE)[1], winslash='/')
code_dir <- normalizePath('.',winslash='/')
repo <- Sys.getenv('MORANK_REPO_DIR', dirname(code_dir))
parent <- normalizePath(Sys.getenv('MORANK_INPUT_DIR'),winslash='/',mustWork=TRUE)
result_dir <- file.path(project, 'results/final_model')
dir.create(result_dir, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table); library(Matrix); library(xgboost)})
setDTthreads(1L)
input <- function(path) file.path(parent,path)
read_result <- function(name) fread(file.path(result_dir,name))
save_result <- function(data,name) fwrite(data,file.path(result_dir,name),bom=TRUE)
save_object <- function(data,name) saveRDS(data,file.path(result_dir,name))
load_object <- function(name) readRDS(file.path(result_dir,name))
read_input <- function(name) fread(input(name))
finite_bh <- function(p) p.adjust(ifelse(is.finite(p),p,1),'BH')
with_directory <- function(path,expr) {
  previous <- setwd(path)
  on.exit(setwd(previous))
  force(expr)
}
focus_genes <- c('ABCC3','NQO1','IRS2','KIT','FGFR3')
discovery_ranks <- function() read_input('round33_model_optimization/results/TCGA_exploratory_ranks.csv')
candidate_list <- function() discovery_ranks()[final_candidate==TRUE][order(final_rank),.(gene,final_rank)]
canonical_map <- function() read_input('round9/results/gene_mapping_audit.csv')
candidate_mapping <- function() merge(candidate_list(),canonical_map()[mapping_status %in%
  c('exact_approved','unique_previous'),.(gene=input_gene,canonical_symbol)],by='gene',sort=FALSE)

load_frozen_rule <- function() {
  with_directory(file.path(repo,'simulation/robustness'),{
    source('adaptive.R',local=.GlobalEnv)
    gate <- read_frozen_model('models/pool_gate.json')
    spec <- fread('settings/frozen_specialist_MoRank37.csv')
    specialist <- read_frozen_model(file.path('models',paste0(spec$model,'.json')))
    i <- which(choices$family=='MoRank37')
    old <- old_models[[i]]
    mixed <- new_models[[i]]
    function(counts,features) {
      prob <- predict(gate,xgb.DMatrix(matrix(pool_profile(counts),nrow=1),nthread=1))
      original <- predict_features(old$model,features,old$cols)
      revised <- if(prob>=.5) predict_features(mixed,features,old$cols) else
        combination(predict_features(specialist,features,old$cols),original,spec$weight)
      list(original=original,score=revised,pln_probability=prob)
    }
  })
}

tcga_data <- function() {
  path <- Sys.getenv('MORANK_TCGA_RDATA',input('TCGA-LUAD_list.RData'))
  if(!file.exists(path)) stop('Missing TCGA input: ', path)
  env <- new.env()
  load(path,envir=env)
  counts <- as.matrix(env$list_data$puried_data)
  group <- as.character(env$list_data$samplename_type)
  if(!all(group %in% c('0','1'))) {
    code <- substr(colnames(counts),14,15)
    group <- ifelse(code %in% c('01','02','03'),'1',ifelse(code %in% c('10','11'),'0',NA))
  }
  stopifnot(sum(group=='1')==541L,sum(group=='0')==59L)
  counts[!is.finite(counts)|counts<0] <- 0
  ids <- trimws(rownames(counts))
  keep <- !is.na(ids)&nzchar(ids)
  counts <- counts[keep,,drop=FALSE]
  rownames(counts) <- ids[keep]
  if(anyDuplicated(rownames(counts))) counts <- rowsum(counts,rownames(counts),reorder=FALSE)
  counts <- round(counts)
  counts <- counts[rowSums(counts>=10)>=30L,,drop=FALSE]
  genes <- read_input('round7/results/all_gene_complement_audit.csv')$gene
  stopifnot(nrow(counts)==29196L,setequal(genes,rownames(counts)))
  counts <- counts[genes,,drop=FALSE]
  expr <- log2(sweep(counts,2,pmax(colSums(counts),1),'/')*1e6+1)
  list(counts=counts,expr=expr,group=group)
}
