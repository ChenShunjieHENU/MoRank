suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
args <- commandArgs(TRUE)
root <- normalizePath(args[1], winslash="/")
out <- function(x, name) fwrite(x,file.path(root,"results",paste0(name,".csv")))
read_workers <- function(directory, pattern, n) {
  files <- list.files(file.path(root,directory),pattern,full.names=TRUE)
  stopifnot(length(files)==n)
  rbindlist(lapply(files,fread))
}
stratified <- function(values, key, column) {
  cells <- values[, .(m = mean(get(column)), v = var(get(column)) / .N, n = .N),
                  by = c(key, 'design', 'mechanism')]
  cells[, {
    centre <- mean(m)
    se <- sqrt(sum(v)) / .N
    df <- sum(v)^2 / sum(v^2 / (n - 1))
    radius <- if (se > 0) qt(.975, df) * se else 0
    .(auc = centre, ci_low = centre - radius, ci_high = centre + radius,
      p = if (se > 0) 2 * pt(-abs(centre / se), df) else if (centre == 0) 1 else 0)
  }, by = key]
}

if (args[2]=="stress") {
  replay <- read_workers("null","^replay_[0-9]+[.]csv$",12L)
  stopifnot(nrow(replay)==640L)
  old <- fread(file.path(root,"primary","results","stress_auc.csv"))
  stopifnot(uniqueN(old$seed)==80L)
  old <- old[!method %in% c("WSKR_A","ClassicalStack","OriginalWSKR"),
              .(design,mechanism=match(scenario,c("S1","S2","S3","S4","S6","S7","S8","S9")),rep,seed,n1,n0,method,auc)]
  new <- replay[suite=="Stress",.(design,mechanism,rep,seed,n1,n0,method,auc)]
  stopifnot(setequal(old$seed,new$seed))
  complete <- rbind(new,old,use.names=TRUE)
  methods <- stratified(complete,"method","auc")
  pairs <- merge(complete[!method %in% c("MoRank","Wilcoxon")],complete[method=="MoRank",.(seed,target=auc)],by="seed")
  pairs[,gain:=target-auc]
  contrasts <- stratified(pairs,"method","gain")
  stopifnot(nrow(contrasts)==13L)
  contrasts[,q:=p.adjust(p,"BH")]
  setorder(methods,-auc)
  out(methods,"final_figure_stress_macro")
  out(contrasts,"final_figure_stress_contrasts")
} else {
  stopifnot(args[2]=="null")
  null <- read_workers("null","^calibration_[0-9]+[.]csv$",26L)
  stopifnot(nrow(null) == 104L, all(null$B == 99L), all(null$genes == 2000L),
            uniqueN(null, by = c('seed','method')) == 104L)
  null[, rate := rejected / genes]
  out(null, 'final_figure_null_datasets')
  intervals <- null[, {
    mean_rate <- mean(rate)
    radius <- qt(.975, .N - 1) * sd(rate) / sqrt(.N)
    .(rate = mean_rate, ci_low = max(0, mean_rate - radius), ci_high = min(1, mean_rate + radius),
      n = .N, routed_pln = sum(route_probability >= .5))
  }, by = .(stage, family, design, n1, n0, method)]
  stopifnot(all(intervals[stage == 'Initial', n] == 6L), all(intervals[stage == 'Follow-up', n] == 8L))
  out(intervals, 'final_figure_null_summary')
}
