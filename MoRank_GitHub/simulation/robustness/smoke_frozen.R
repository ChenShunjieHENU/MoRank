source("adaptive.R")
set.seed(9401)
for (family in c("NB", "PLN")) {
  sim <- simulate_case(family, 1L, 1L, 9401L)
  features <- both_features(sim$counts, sim$n1, sim$n0)
  gate <- read_frozen_model("models/pool_gate.json")
  probability <- predict(gate, xgb.DMatrix(matrix(pool_profile(sim$counts), nrow = 1), nthread = 1))
  stopifnot(ncol(features) == 48L, is.finite(probability),
            (probability >= .5) == (family == "PLN"))
  for (i in seq_len(nrow(choices))) {
    pair <- predict_new_old(i, features)
    spec <- fread(file.path("settings", paste0("frozen_specialist_", choices$family[i], ".csv")))
    specialist <- read_frozen_model(file.path("models", paste0(spec$model, ".json")))
    nb <- predict_features(specialist, features, features_for(choices$family[i], colnames(features)))
    score <- if (probability >= .5) pair$new else combination(nb, pair$old, spec$weight)
    stopifnot(length(score) == nrow(sim$counts), all(is.finite(score)))
    cat(family, choices$family[i], "AUC", auc(score, sim$truth), "\n")
  }
}
