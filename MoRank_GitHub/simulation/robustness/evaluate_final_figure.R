# Execute from the frozen robustness benchmark directory.
source("adaptive.R")

arguments <- commandArgs(trailingOnly = TRUE)
phase <- arguments[1]
worker <- as.integer(arguments[2])
workers <- as.integer(arguments[3])
output <- normalizePath(arguments[4], mustWork = TRUE)
stopifnot(phase %in% c("replay", "calibration"), worker %in% seq_len(workers))

gate <- read_frozen_model("models/pool_gate.json")
specialists <- lapply(choices$family, function(family) {
  row <- fread(file.path("settings", paste0("frozen_specialist_", family, ".csv")))
  list(model = read_frozen_model(file.path("models", paste0(row$model, ".json"))),
       weight = row$weight)
})
route <- function(counts) {
  predict(gate, xgb.DMatrix(matrix(pool_profile(counts), nrow = 1), nthread = 1))
}
normalize_counts <- function(counts) {
  log2(sweep(counts, 2, pmax(colSums(counts), 1), "/") * 1e6 + 1)
}
score <- function(features, pln_prob) {
  values <- lapply(seq_len(nrow(choices)), function(index) {
    columns <- features_for(choices$family[index], colnames(features))
    if (pln_prob >= .5) {
      return(predict_features(new_models[[index]], features, columns))
    }
    branch <- specialists[[index]]
    base <- old_models[[index]]
    combination(predict_features(branch$model, features, columns),
                predict_features(base$model, features, base$cols), branch$weight)
  })
  matrix <- do.call(cbind, values)
  colnames(matrix) <- ifelse(choices$family == "Classical", "Classical", "MoRank")
  matrix
}

if (phase == "replay") {
  jobs <- rbindlist(lapply(c(Primary = "final", Stress = "stress"), function(name) {
    files <- list.files(file.path(Sys.getenv("MORANK_PRIMARY_DIR"), paste0(name, "_raw")), "\\.rds$", full.names = TRUE)
    fields <- strcapture("^d([0-9]+)_s([0-9]+)_r([0-9]+)\\.rds$", basename(files),
                         list(design = integer(), mechanism = integer(), rep = integer()))
    data.table(suite = if (name == "final") "Primary" else "Stress", file = files, fields)
  }))
  stopifnot(nrow(jobs) == 320L, !anyNA(jobs))
  result <- lapply(seq(worker, nrow(jobs), by = workers), function(index) {
    job <- jobs[index]
    saved <- readRDS(job$file)
    family <- if (job$suite == "Primary") "NB" else "Sensitivity"
    sim <- simulate_case(family, job$design, job$mechanism, saved$seed)
    stopifnot(identical(as.logical(sim$truth), as.logical(saved$truth)))
    features <- distribution_features(normalize_counts(sim$counts),
                                      sim$n0 + seq_len(sim$n1), seq_len(sim$n0))
    reference <- auc(saved$scores[, "KS"], sim$truth)
    stopifnot(abs(auc(features[, "KS"], sim$truth) - reference) < 1e-9)
    p <- route(sim$counts)
    predicted <- score(features, p)
    data.table(suite = job$suite, design = job$design, mechanism = job$mechanism,
               rep = job$rep, seed = saved$seed, n1 = sim$n1, n0 = sim$n0, route_probability = p,
               method = colnames(predicted), auc = apply(predicted, 2, auc, truth = sim$truth))
  })
} else {
  initial <- as.data.table(expand.grid(family = c("NB", "PLN"), design = 1:3, replicate = 1:6))
  followup <- as.data.table(expand.grid(family = c("NB", "PLN"), design = 3L, replicate = 1:8))
  jobs <- rbindlist(list(initial[, stage := "Initial"], followup[, stage := "Follow-up"]))
  jobs[, seed := 882000000L + ifelse(stage == "Initial", 100000L, 200000L) +
         match(family, c("NB", "PLN")) * 10000L + design * 100L + replicate]
  stopifnot(nrow(jobs) == 52L, uniqueN(jobs$seed) == 52L)
  result <- lapply(seq(worker, nrow(jobs), by = workers), function(index) {
    job <- jobs[index]
    sizes <- primary_designs[job$design, ]
    n1 <- sizes[1]
    n0 <- sizes[2]
    if (job$family == "PLN") {
      counts <- generate_pln(1L, n1, n0, job$seed, signal_n = 0L)$counts
    } else {
      set.seed(job$seed)
      config <- default_config()
      config$single_signal_n <- 0L
      sim <- simulate_pool_fast("S1", max(n1, n0), config)
      stopifnot(!any(sim$truth$is_signal))
      counts <- sim$counts[, c(seq_len(n0), max(n1, n0) + seq_len(n1)), drop = FALSE]
    }
    x <- normalize_counts(counts)
    p <- route(counts)
    n <- n0 + n1
    evaluate <- function(case) score(distribution_features(x, case, setdiff(seq_len(n), case)), p)
    observed <- evaluate(n0 + seq_len(n1))
    B <- 99L
    exceed <- matrix(0L, nrow(x), ncol(observed))
    set.seed(job$seed + 1000000L)
    for (permutation in seq_len(B)) {
      case <- sort(sample.int(n, n1))
      exceed <- exceed + (evaluate(case) >= observed)
    }
    permutation_p <- (1 + exceed) / (B + 1)
    cat(job$stage, job$family, job$design, job$seed, "complete\n")
    flush.console()
    data.table(stage = job$stage, family = job$family, design = job$design,
               n1 = n1, n0 = n0, seed = job$seed, route_probability = p, B = B,
               genes = nrow(x), method = colnames(observed), rejected = colSums(permutation_p <= .05))
  })
}
fwrite(rbindlist(result), file.path(output, sprintf("%s_%02d.csv", phase, worker)))
