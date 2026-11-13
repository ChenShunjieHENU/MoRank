source_original <- function(path) {
  previous <- getwd()
  on.exit(setwd(previous))
  setwd(path)
  source("original_simulation.R", local = .GlobalEnv)
  source("simulate_pool_fast.R", local = .GlobalEnv)
  source("redesign_core.R", local = .GlobalEnv)
  source("model_io.R", local = .GlobalEnv)
}
source_original("../reference")
source("cross_generator_core.R")

local_features <- function(x, case, control) {
  n <- ncol(x)
  ne <- length(case) * length(control) / n
  ranks <- (matrixStats::rowRanks(x, ties.method="average", useNames=FALSE)-.5)/n
  groups <- lapply(c(.1,.2,.4), function(width) {
    endpoints <- seq(0,1-width,by=width/2)
    values <- vapply(endpoints, function(start) {
      block <- (ranks > start) & (ranks <= start + width + 1e-12)
      a <- rowMeans(block[,case,drop=FALSE])
      b <- rowMeans(block[,control,drop=FALSE])
      p <- rowMeans(block)
      positive_ratio(ne*(a-b)^2, n*p*(1-p)/(n-1))
    }, numeric(nrow(x)))
    cbind(matrixStats::rowMaxs(values), rowMeans(values))
  })
  windows <- do.call(cbind, groups)
  colnames(windows) <- paste0("Interval",rep(c(10,20,40),each=2),c("Max","Mean"))
  p <- c(.05,.15,.25,.4,.5,.6,.75,.85,.95)
  q1 <- matrixStats::rowQuantiles(x[,case,drop=FALSE],probs=p,type=8,useNames=FALSE)
  q0 <- matrixStats::rowQuantiles(x[,control,drop=FALSE],probs=p,type=8,useNames=FALSE)
  spread <- function(q) cbind(q[,9]-q[,1], q[,7]-q[,3],
                             q[,3]-q[,1], q[,9]-q[,7])
  scale <- matrixStats::rowSds(x)
  floor <- pmax(0.02*scale,1e-8)
  spreads <- abs(log((spread(q1)+floor)/(spread(q0)+floor)))
  colnames(spreads) <- c("Central90Ratio","IQRRatio","LowerTailRatio","UpperTailRatio")
  delta <- q1-q0
  delta <- delta-rowMeans(delta)
  shape <- positive_ratio(ne*rowMeans(delta^2), scale^2)
  result <- log1p(cbind(windows,spreads,QuantileContrastShape=shape))
  stopifnot(all(is.finite(result)), !anyDuplicated(colnames(result)))
  result
}

both_features <- function(counts, n1, n0) {
  x <- log2(sweep(counts,2,pmax(colSums(counts),1),"/")*1e6+1)
  case <- n0+seq_len(n1)
  control <- seq_len(n0)
  cbind(distribution_features(x,case,control),local_features(x,case,control))
}

mechanisms <- c("S1","S2","S3","S4","S6","S7","S8","S9")
primary_designs <- rbind(c(50L,50L),c(200L,200L),c(541L,59L))
small_designs <- rbind(c(20L,20L),c(80L,80L),c(120L,30L))
cells_for <- function(suite) {
  if (suite=="NB") expand.grid(design=1:3, mechanism=1:8)
  else if (suite=="Sensitivity") expand.grid(design=1:2, mechanism=1:8)
  else expand.grid(design=1:3, mechanism=1:4)
}
sample_sizes <- function(suite,design) {
  sizes <- switch(suite,NB=primary_designs,Sensitivity=primary_designs[c(1,3),],
                  PLN=primary_designs,PLN_small=small_designs)
  sizes[design,]
}
stress_config <- function(varied) {
  c <- default_config()
  c$tail_pi <- .02
  c$tail_ratio <- 8
  c$dispersion_ratio <- 1.4
  c$skew_sigma_case <- .4
  c$proportion_case <- .5
  c$modality_low_ratio <- .75
  c$modality_high_ratio <- 1.25
  c$mean_log2fc_min <- .3
  c$mean_log2fc_max <- .8
  if (varied) {
    c$tail_pi <- runif(1,.015,.035)
    c$tail_ratio <- runif(1,6,9)
    c$dispersion_ratio <- runif(1,1.25,1.65)
    c$skew_sigma_case <- runif(1,.32,.52)
    c$proportion_case <- runif(1,.42,.58)
    width <- runif(1,.17,.34)
    c$modality_low_ratio <- 1-width
    c$modality_high_ratio <- 1+width
  }
  c
}
simulate_case <- function(suite, design, mechanism, seed, varied=FALSE) {
  n <- sample_sizes(suite,design)
  n1 <- n[[1]]
  n0 <- n[[2]]
  if (startsWith(suite,"PLN")) return(generate_pln(mechanism,n1,n0,seed))
  set.seed(seed)
  config <- if(suite=="Sensitivity") stress_config(varied) else default_config()
  max_n <- max(n)
  sim <- simulate_pool_fast(mechanisms[mechanism],max_n,config)
  cols <- c(seq_len(n0),max_n+seq_len(n1))
  list(counts=sim$counts[,cols,drop=FALSE],truth=sim$truth$is_signal,n1=n1,n0=n0)
}

suites <- c("NB","Sensitivity","PLN","PLN_small")
jobs_for <- function(split) {
  replicates <- c(fit=4L,tune=2L,test=3L,test2=3L)[[split]]
  base <- c(fit=810000000L,tune=820000000L,test=830000000L,test2=840000000L)[[split]]
  data.table::rbindlist(lapply(seq_along(suites),function(i) {
    z <- data.table::as.data.table(cells_for(suites[i]))
    z <- z[rep(seq_len(.N),each=replicates)]
    z[,replicate:=rep(seq_len(replicates),times=.N/replicates)]
    z[,suite:=suites[i]]
    z[,seed:=base+i*100000L+design*1000L+mechanism*10L+replicate]
    z
  }))
}
