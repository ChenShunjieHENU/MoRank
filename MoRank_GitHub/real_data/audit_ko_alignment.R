suppressPackageStartupMessages({library(data.table);library(Matrix);library(RSpectra);library(cli)})
args <- commandArgs(TRUE)
root <- args[1]; worker <- as.integer(args[2]); workers <- as.integer(args[3])
src <- file.path(Sys.getenv("MORANK_INPUT_DIR"),"round10")
ref <- readRDS(file.path(src,"data/IRS2_reference_subnetwork.rds"))
N <- ref$network
old <- readRDS(file.path(src,"results/reference_KO/perturbations_all102.rds"))
targets <- c("IRS2","SHAM",ref$controls)
stopifnot(length(targets)==102,setequal(targets,names(old)))
algorithm <- Sys.getenv("MORANK_ALIGNMENT_R")
stopifnot(file.exists(algorithm))
align <- eval(parse(file=algorithm))
dest <- file.path(root,"common/results",sprintf("ko_align_%02d.csv",worker))
rows <- list()
for (k in which((seq_along(targets)-1) %% workers + 1 == worker)) {
 g <- targets[k]; K <- N
 if (g!="SHAM") K[g,] <- 0
 set.seed(1)
 A <- align(N,K,d=2,nCores=1)
 stopifnot(all(dim(A)==c(4000,2)),all(is.finite(A)))
 d <- setNames(sqrt(rowSums((A[1:2000,,drop=FALSE]-A[2001:4000,,drop=FALSE])^2)),ref$nodes)
 prev <- old[[g]]$distances[ref$nodes]
 endpoints <- c(global=sum(d[ref$response_genes]),vapply(ref$modules,function(m)mean(d[m]),numeric(1)))
 prev_endpoints <- c(global=sum(prev[ref$response_genes]),vapply(ref$modules,function(m)mean(prev[m]),numeric(1)))
 rows[[length(rows)+1L]] <- cbind(data.table(target=g,max_abs_distance_error=max(abs(d-prev))),
  as.data.table(as.list(endpoints)),as.data.table(as.list(setNames(prev_endpoints,paste0("old_",names(prev_endpoints))))))
 fwrite(rbindlist(rows),dest)
 saveRDS(d,file.path(root,"common/results",sprintf("ko_distance_%03d.rds",k)))
 cat("KO",worker,k,g,signif(max(abs(d-prev)),4),"\n"); flush.console()
}
