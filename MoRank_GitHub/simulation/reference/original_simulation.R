
# ============================================================
# WSKR Simulation
# Negative-binomial simulation and AUC evaluation
# ============================================================

library(data.table)

# 1. 参数 ---------------------------------------------------

default_config <- function() {
  list(
    out = "results",
    p = 2000L,
    reps = 30L,
    sample_sizes = c(50L, 100L, 200L),
    
    calibration_reps = 10L,
    calibration_n = 100L,
    alpha_grid = c(0, 0.25, 0.50, 0.75, 1),
    alpha_fixed = NA_real_,
    
    seed = 20260916L,
    calibration_seed = 20260917L,
    
    q_log_mean = log(20),
    q_log_sd = 1.25,
    q_min = 0.01,
    q_max = 5000,
    
    phi_log_mean = log(0.08),
    phi_log_sd = 0.80,
    phi_min = 1e-4,
    phi_max = 3,
    
    library_log_sd = 0.45,
    
    mean_log2fc_min = 0.60,
    mean_log2fc_max = 1.50,
    
    skew_sigma_control = 0.10,
    skew_sigma_case = 0.62,
    
    tail_pi = 0.045,
    tail_ratio = 5.5,
    tail_dispersion = 3,
    
    dispersion_ratio = 2.0,
    
    proportion_control = 0.25,
    proportion_case = 0.65,
    proportion_high_ratio = 1.40,
    
    modality_low_ratio = 0.60,
    modality_high_ratio = 1.40,
    
    single_signal_n = 120L,
    mixed_signal_n = 40L,
    score_eps = 1e-8,
    
    skip_count_methods = TRUE
  )
}

scenario_levels <- c(
  "S0", "S1", "S2", "S3", "S4",
  "S6", "S7", "S8", "S9"
)

scenario_labels <- c(
  S0 = "Null",
  S1 = "Location",
  S2 = "Skewness",
  S3 = "Tail",
  S4 = "Compound",
  S6 = "Dispersion",
  S7 = "Proportion",
  S8 = "Modality",
  S9 = "Mechanism-diverse mixture"
)

formal_labels <- c(
  S0 = "Null reference scenario",
  S1 = "Location-shift scenario",
  S2 = "Skewness-altered scenario",
  S3 = "Tail-altered scenario",
  S4 = "Compound shape-altered scenario",
  S6 = "Dispersion-altered scenario",
  S7 = "Mixture-proportion-altered scenario",
  S8 = "Modality-altered scenario",
  S9 = "Mechanism-diverse mixture scenario"
)

method_levels <- c(
  "WSKR", "SR", "KR", "KS",
  "WaddR_shape", "Wilcoxon", "DESeq2", "edgeR"
)

scenario_name <- function(x) unname(scenario_labels[x])


# 2. 负二项数据生成 -----------------------------------------

clip_log_normal <- function(n, meanlog, sdlog, lower, upper) {
  x <- exp(rnorm(n, meanlog, sdlog))
  pmin(pmax(x, lower), upper)
}

normalize_mean <- function(x) {
  x / max(mean(x), 1e-12)
}

simulate_library <- function(n, config) {
  s <- exp(rnorm(
    n,
    -0.5 * config$library_log_sd^2,
    config$library_log_sd
  ))
  s / mean(s)
}

nb_draw <- function(mu, phi) {
  rnbinom(
    length(mu),
    mu = pmax(mu, 1e-8),
    size = 1 / pmax(phi, 1e-8)
  )
}

tail_draw <- function(mu, phi, config) {
  pi_tail <- config$tail_pi
  high <- rbinom(length(mu), 1, pi_tail) == 1
  
  mu_high <- config$tail_ratio * mu
  mu_low <- (mu - pi_tail * mu_high) / (1 - pi_tail)
  
  mu_mix <- ifelse(high, mu_high, mu_low)
  phi_mix <- ifelse(
    high,
    phi * config$tail_dispersion,
    phi
  )
  
  nb_draw(mu_mix, phi_mix)
}

proportion_draw <- function(mu, phi, pi_high, config, n) {
  high <- rbinom(n, 1, pi_high) == 1
  
  r <- config$proportion_high_ratio
  low <- (1 - pi_high * r) / (1 - pi_high)
  
  mu_i <- ifelse(high, r * mu, low * mu)
  nb_draw(mu_i, phi)
}

modality_draw <- function(mu, phi, config, n) {
  high <- rbinom(n, 1, 0.5) == 1
  
  mu_i <- ifelse(
    high,
    config$modality_high_ratio * mu,
    config$modality_low_ratio * mu
  )
  
  nb_draw(mu_i, phi)
}

generate_parameters <- function(config) {
  data.table(
    gene = sprintf("gene_%05d", seq_len(config$p)),
    base_mean = clip_log_normal(
      config$p,
      config$q_log_mean,
      config$q_log_sd,
      config$q_min,
      config$q_max
    ),
    dispersion = clip_log_normal(
      config$p,
      config$phi_log_mean,
      config$phi_log_sd,
      config$phi_min,
      config$phi_max
    )
  )
}


# 3. 定义真实信号基因 ---------------------------------------

make_truth <- function(scenario, config) {
  
  truth <- data.table(
    gene = sprintf("gene_%05d", seq_len(config$p)),
    signal_type = "Null"
  )
  
  single_types <- c(
    S1 = "Location",
    S2 = "Skewness",
    S3 = "Tail",
    S4 = "Compound",
    S6 = "Dispersion",
    S7 = "Proportion",
    S8 = "Modality"
  )
  
  if (scenario %in% names(single_types)) {
    n <- min(config$single_signal_n, config$p)
    
    truth$signal_type[seq_len(n)] <-
      single_types[[scenario]]
  }
  
  if (scenario == "S9") {
    types <- c(
      "Location", "Skewness", "Tail",
      "Dispersion", "Proportion", "Modality"
    )
    
    n <- min(
      config$mixed_signal_n,
      floor(config$p / length(types))
    )
    
    idx <- 1L
    
    for (type in types) {
      truth$signal_type[idx:(idx + n - 1L)] <- type
      idx <- idx + n
    }
  }
  
  truth[, is_signal := signal_type != "Null"]
  
  truth[, is_nonlocation := signal_type %in% c(
    "Skewness", "Tail", "Compound",
    "Dispersion", "Proportion", "Modality"
  )]
  
  truth
}


# 4. 按信号机制模拟单个基因 ---------------------------------

simulate_one_gene <- function(
    type, q, phi, lib_control, lib_case, config) {
  
  control_mu <- lib_control * q
  case_mu <- lib_case * q
  n <- length(lib_control)
  
  # 位置偏移
  if (type == "Location") {
    beta <- runif(
      1,
      config$mean_log2fc_min,
      config$mean_log2fc_max
    ) * sample(c(-1, 1), 1)
    
    case_mu <- case_mu * 2^beta
    
    return(c(
      nb_draw(control_mu, phi),
      nb_draw(case_mu, phi)
    ))
  }
  
  # 偏度和复合形态改变
  if (type %in% c("Skewness", "Compound")) {
    control_mu <- control_mu * normalize_mean(
      exp(rnorm(
        n,
        -0.5 * config$skew_sigma_control^2,
        config$skew_sigma_control
      ))
    )
    
    case_mu <- case_mu * normalize_mean(
      exp(rnorm(
        n,
        -0.5 * config$skew_sigma_case^2,
        config$skew_sigma_case
      ))
    )
  }
  
  control <- if (type == "Proportion") {
    proportion_draw(
      control_mu, phi,
      config$proportion_control, config, n
    )
  } else {
    nb_draw(control_mu, phi)
  }
  
  case <- switch(
    type,
    Tail = tail_draw(case_mu, phi, config),
    Compound = tail_draw(case_mu, phi, config),
    Dispersion = nb_draw(
      case_mu, phi * config$dispersion_ratio
    ),
    Proportion = proportion_draw(
      case_mu, phi,
      config$proportion_case, config, n
    ),
    Modality = modality_draw(
      case_mu, phi, config, n
    ),
    nb_draw(case_mu, phi)
  )
  
  c(control, case)
}


# 5. 生成完整仿真数据集 -------------------------------------

simulate_dataset <- function(scenario, n_per_group, config) {
  
  gene_params <- generate_parameters(config)
  truth <- make_truth(scenario, config)
  
  group <- c(
    rep("control", n_per_group),
    rep("case", n_per_group)
  )
  
  lib <- simulate_library(2L * n_per_group, config)
  
  control_lib <- lib[seq_len(n_per_group)]
  case_lib <- lib[
    (n_per_group + 1L):(2L * n_per_group)
  ]
  
  counts <- matrix(
    0L,
    nrow = config$p,
    ncol = 2L * n_per_group,
    dimnames = list(
      gene_params$gene,
      paste0("sample_", seq_len(2L * n_per_group))
    )
  )
  
  for (g in seq_len(config$p)) {
    counts[g, ] <- simulate_one_gene(
      truth$signal_type[g],
      gene_params$base_mean[g],
      gene_params$dispersion[g],
      control_lib,
      case_lib,
      config
    )
  }
  
  list(
    counts = counts,
    group = group,
    truth = truth,
    gene_params = gene_params
  )
}


# 6. 计算偏度和峰度 -----------------------------------------

row_skew <- function(x) {
  m <- rowMeans(x)
  z <- sweep(x, 1, m, "-")
  
  s <- sqrt(rowMeans(z^2))
  s[s < 1e-8] <- NA_real_
  
  out <- rowMeans(sweep(z, 1, s, "/")^3)
  out[!is.finite(out)] <- 0
  
  out
}

row_tail <- function(x) {
  m <- rowMeans(x)
  z <- sweep(x, 1, m, "-")
  
  s <- sqrt(rowMeans(z^2))
  s[s < 1e-8] <- NA_real_
  
  out <- rowMeans(sweep(z, 1, s, "/")^4) - 3
  out[!is.finite(out)] <- 0
  
  sign(out) * log1p(abs(out))
}


# 7. Wasserstein shape --------------------------------------

wasserstein_shape <- function(
    control, case,
    probs = seq(0.01, 0.99, length.out = 99L)) {
  
  qc <- as.numeric(quantile(
    control, probs = probs, names = FALSE, type = 8
  ))
  
  qt <- as.numeric(quantile(
    case, probs = probs, names = FALSE, type = 8
  ))
  
  sd_c <- sd(control)
  sd_t <- sd(case)
  
  rho <- suppressWarnings(cor(qc, qt))
  
  if (!is.finite(rho) ||
      !is.finite(sd_c) ||
      !is.finite(sd_t)) return(0)
  
  max(0, 2 * sd_c * sd_t * (1 - rho))
}


# 8. 各方法的基因评分 ---------------------------------------

compute_scores <- function(counts, group, alpha) {
  
  lib <- colSums(counts)
  lib[lib <= 0] <- 1
  
  log_expr <- log2(
    t(t(counts) / lib) * 1e6 + 1
  )
  
  control <- log_expr[
    , group == "control", drop = FALSE
  ]
  
  case <- log_expr[
    , group == "case", drop = FALSE
  ]
  
  skew_c <- row_skew(control)
  skew_t <- row_skew(case)
  
  tail_c <- row_tail(control)
  tail_t <- row_tail(case)
  
  dskew <- skew_t - skew_c
  dtail <- tail_t - tail_c
  
  sr <- dskew / (
    abs(skew_t) + abs(skew_c) + 1e-8
  )
  
  kr <- dtail / (
    abs(tail_t) + abs(tail_c) + 1e-8
  )
  
  wskr <- (
    alpha * dskew + (1 - alpha) * dtail
  ) / (
    alpha * (abs(skew_t) + abs(skew_c)) +
      (1 - alpha) * (abs(tail_t) + abs(tail_c)) +
      1e-8
  )
  
  ks <- vapply(
    seq_len(nrow(log_expr)),
    function(i) {
      suppressWarnings(
        ks.test(control[i, ], case[i, ])$statistic[[1]]
      )
    },
    numeric(1)
  )
  
  wshape <- vapply(
    seq_len(nrow(log_expr)),
    function(i) {
      wasserstein_shape(control[i, ], case[i, ])
    },
    numeric(1)
  )
  
  wilcox <- vapply(
    seq_len(nrow(log_expr)),
    function(i) {
      suppressWarnings(
        wilcox.test(
          control[i, ],
          case[i, ],
          exact = FALSE
        )$p.value
      )
    },
    numeric(1)
  )
  
  scores <- data.table(
    gene = rownames(counts),
    
    delta_skew = dskew,
    delta_tail = dtail,
    
    skew_case = skew_t,
    skew_control = skew_c,
    
    tail_case = tail_t,
    tail_control = tail_c,
    
    WSKR = abs(wskr),
    SR = abs(sr),
    KR = abs(kr),
    KS = ks,
    WaddR_shape = wshape,
    
    Wilcoxon = -log10(
      pmax(wilcox, .Machine$double.xmin)
    )
  )
  
  methods <- c(
    "WSKR", "SR", "KR", "KS",
    "WaddR_shape", "Wilcoxon"
  )
  
  for (method in methods) {
    set(
      scores,
      which(!is.finite(scores[[method]])),
      method,
      0
    )
  }
  
  list(scores = scores, log_expr = log_expr)
}


# 9. DESeq2 和 edgeR（可选） --------------------------------

run_deseq2 <- function(counts, group) {
  
  coldata <- data.frame(
    condition = factor(
      group, levels = c("control", "case")
    ),
    row.names = colnames(counts)
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    round(counts),
    coldata,
    design = ~ condition
  )
  
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  res <- as.data.frame(DESeq2::results(
    dds,
    contrast = c("condition", "case", "control")
  ))
  
  out <- data.table(
    gene = rownames(res),
    DESeq2 = -log10(
      pmax(res$padj, .Machine$double.xmin)
    )
  )
  
  out[!is.finite(DESeq2), DESeq2 := 0]
  out
}

run_edger <- function(counts, group) {
  
  dge <- edgeR::DGEList(
    counts = round(counts),
    group = factor(group)
  )
  
  keep <- edgeR::filterByExpr(
    dge,
    group = dge$samples$group
  )
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ dge$samples$group)
  
  dge <- edgeR::estimateDisp(dge, design)
  
  fit <- edgeR::glmQLFit(dge, design)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)
  
  tt <- edgeR::topTags(
    qlf,
    n = Inf,
    sort.by = "none"
  )$table
  
  out <- data.table(
    gene = rownames(tt),
    edgeR = -log10(
      pmax(tt$FDR, .Machine$double.xmin)
    )
  )
  
  out[!is.finite(edgeR), edgeR := 0]
  out
}


# 10. AUC 计算 ----------------------------------------------

auc_score <- function(score, truth) {
  
  ok <- is.finite(score) & !is.na(truth)
  
  score <- score[ok]
  truth <- as.logical(truth[ok])
  
  p <- sum(truth)
  n <- sum(!truth)
  
  if (p == 0 || n == 0) return(NA_real_)
  
  r <- rank(score, ties.method = "average")
  
  (sum(r[truth]) - p * (p + 1) / 2) / (p * n)
}


# 11. 单次仿真评价 ------------------------------------------

evaluate_iteration <- function(
    scenario, n_per_group, rep_id, config, alpha) {
  
  sim <- simulate_dataset(
    scenario, n_per_group, config
  )
  
  scored <- compute_scores(
    sim$counts, sim$group, alpha
  )
  
  scores <- scored$scores
  
  if (!config$skip_count_methods) {
    
    scores <- merge(
      scores,
      run_deseq2(sim$counts, sim$group),
      by = "gene", all.x = TRUE
    )
    
    scores <- merge(
      scores,
      run_edger(sim$counts, sim$group),
      by = "gene", all.x = TRUE
    )
  }
  
  scores <- merge(
    scores,
    sim$truth[, .(
      gene, signal_type, is_signal, is_nonlocation
    )],
    by = "gene",
    sort = FALSE
  )
  
  target_list <- list(
    Any = scores$is_signal,
    Nonlocation = scores$is_nonlocation,
    Location = scores$signal_type == "Location",
    Skewness = scores$signal_type == "Skewness",
    Tail = scores$signal_type == "Tail",
    Compound = scores$signal_type == "Compound",
    Dispersion = scores$signal_type == "Dispersion",
    Proportion = scores$signal_type == "Proportion",
    Modality = scores$signal_type == "Modality",
    Mixed_nonlocation = scores$signal_type %in% c(
      "Skewness", "Tail", "Dispersion",
      "Proportion", "Modality"
    )
  )
  
  methods <- intersect(method_levels, names(scores))
  rows <- list()
  
  for (method in methods) {
    for (target in names(target_list)) {
      
      truth <- target_list[[target]]
      
      if (sum(truth) == 0 || sum(!truth) == 0) next
      
      rows[[length(rows) + 1L]] <- data.table(
        scenario = scenario,
        scenario_name = scenario_name(scenario),
        scenario_formal_name = formal_labels[[scenario]],
        
        n_per_group = n_per_group,
        rep = rep_id,
        
        method = method,
        target = target,
        
        auc = auc_score(
          scores[[method]], truth
        ),
        
        alpha = alpha,
        signal_n = sum(truth),
        total_genes = nrow(scores)
      )
    }
  }
  
  # 无差异场景参考分布
  null_scores <- if (scenario == "S0") {
    
    rbindlist(lapply(methods, function(method) {
      data.table(
        scenario = scenario,
        scenario_name = scenario_name(scenario),
        n_per_group = n_per_group,
        rep = rep_id,
        method = method,
        
        score_q95 = as.numeric(quantile(
          scores[[method]], 0.95, na.rm = TRUE
        )),
        
        score_max = max(
          scores[[method]], na.rm = TRUE
        )
      )
    }))
    
  } else {
    data.table()
  }
  
  list(
    metrics = rbindlist(rows, fill = TRUE),
    null = null_scores
  )
}


# 12. 独立 alpha 校准 ---------------------------------------

run_alpha_calibration <- function(config) {
  
  set.seed(config$calibration_seed)
  
  rows <- list()
  
  for (rep_id in seq_len(config$calibration_reps)) {
    
    for (scenario in c("S2", "S3", "S4")) {
      
      sim <- simulate_dataset(
        scenario,
        config$calibration_n,
        config
      )
      
      scored <- compute_scores(
        sim$counts, sim$group, alpha = 0.5
      )
      
      s <- scored$scores
      truth <- sim$truth$is_signal
      
      for (alpha in config$alpha_grid) {
        
        score <- abs(
          (
            alpha * s$delta_skew +
              (1 - alpha) * s$delta_tail
          ) / (
            alpha * (
              abs(s$skew_case) +
                abs(s$skew_control)
            ) +
              (1 - alpha) * (
                abs(s$tail_case) +
                  abs(s$tail_control)
              ) +
              config$score_eps
          )
        )
        
        rows[[length(rows) + 1L]] <- data.table(
          rep = rep_id,
          scenario = scenario,
          alpha = alpha,
          auc = auc_score(score, truth)
        )
      }
    }
  }
  
  dt <- rbindlist(rows)
  
  summary <- dt[, .(
    auc_mean = mean(auc),
    auc_sd = sd(auc)
  ), by = .(scenario, alpha)]
  
  choice <- summary[
    , .(objective = mean(auc_mean)),
    by = alpha
  ][order(-objective)][1]
  
  list(
    long = dt,
    summary = summary,
    choice = choice
  )
}


# 13. 主程序 ------------------------------------------------

main <- function() {
  
  config <- default_config()
  
  out_dir <- file.path(config$out, "tables")
  
  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # 独立校准
  calibration <- run_alpha_calibration(config)
  
  selected_alpha <- if (is.finite(config$alpha_fixed)) {
    config$alpha_fixed
  } else {
    calibration$choice$alpha
  }
  
  fwrite(
    calibration$long,
    file.path(out_dir, "alpha_calibration_long.csv")
  )
  
  fwrite(
    calibration$summary,
    file.path(out_dir, "alpha_calibration_summary.csv")
  )
  
  fwrite(
    calibration$choice,
    file.path(out_dir, "alpha_selected.csv")
  )
  
  cat("Selected alpha:", selected_alpha, "\n")
  
  # 主仿真
  set.seed(config$seed)
  
  rows <- list()
  null_rows <- list()
  
  for (n_per_group in config$sample_sizes) {
    
    for (scenario in scenario_levels) {
      
      for (rep_id in seq_len(config$reps)) {
        
        result <- evaluate_iteration(
          scenario,
          n_per_group,
          rep_id,
          config,
          selected_alpha
        )
        
        if (nrow(result$metrics)) {
          rows[[length(rows) + 1L]] <- result$metrics
        }
        
        if (nrow(result$null)) {
          null_rows[[length(null_rows) + 1L]] <- result$null
        }
      }
      
      cat(
        "Completed:", scenario,
        "| n =", n_per_group, "\n"
      )
    }
  }
  
  metrics <- rbindlist(rows, fill = TRUE)
  null <- rbindlist(null_rows, fill = TRUE)
  
  # AUC 原始结果
  fwrite(
    metrics,
    file.path(out_dir, "auc_metrics_long.csv")
  )
  
  # AUC 汇总
  auc_summary <- metrics[, .(
    auc_mean = mean(auc),
    auc_sd = sd(auc),
    n_rep = .N
  ), by = .(
    scenario, scenario_name,
    n_per_group, method, target
  )]
  
  fwrite(
    auc_summary,
    file.path(out_dir, "auc_metrics_summary.csv")
  )
  
  # 无差异参考结果
  fwrite(
    null,
    file.path(out_dir, "null_reference_long.csv")
  )
  
  null_summary <- null[, .(
    score_q95_mean = mean(score_q95),
    score_max_mean = mean(score_max),
    n_rep = .N
  ), by = .(
    scenario, scenario_name,
    n_per_group, method
  )]
  
  fwrite(
    null_summary,
    file.path(out_dir, "null_reference_summary.csv")
  )
  
  # 运行参数
  fwrite(
    data.table(
      parameter = names(config),
      value = vapply(
        config, paste, character(1), collapse = ","
      )
    ),
    file.path(out_dir, "simulation_parameters.csv")
  )
  
  cat("\nSimulation completed.\n")
  cat("Results:", out_dir, "\n")
}

if (sys.nframe() == 0L) main()
