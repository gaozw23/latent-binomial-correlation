.phase3_project_root <- function() {
  normalizePath(getOption("latent_binomial_project_root", "."), winslash = "/", mustWork = TRUE)
}

open_phase3_cluster <- function(workers) {
  workers <- as.integer(workers)
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("workers must be a positive integer", call. = FALSE)
  }
  if (workers == 1L) return(NULL)
  root <- .phase3_project_root()
  libraries <- .libPaths()
  cluster <- parallel::makePSOCKcluster(workers, rscript_args = "--no-init-file")
  parallel::clusterExport(cluster, c("root", "libraries"), envir = environment())
  setup <- parallel::clusterEvalQ(cluster, {
    .libPaths(libraries)
    options(latent_binomial_project_root = root)
    files <- sort(list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE))
    invisible(lapply(files, sys.source, envir = .GlobalEnv))
    TRUE
  })
  if (!all(unlist(setup))) {
    parallel::stopCluster(cluster)
    stop("Pairwise-simulation worker initialisation failed", call. = FALSE)
  }
  cluster
}

phase3_relative_path <- function(path, root = .phase3_project_root()) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(normalizePath(root, winslash = "/", mustWork = TRUE), "/")
  sub(paste0("^", gsub("([.()+*?^$|{}\\[\\]\\\\])", "\\\\\\1", prefix)), "", path)
}

generate_phase3_stream_files <- function(root = .phase3_project_root()) {
  scenarios <- read_pairwise_scenarios(file.path(root, "config", "pairwise_scenarios.csv"))
  full_map <- do.call(rbind, lapply(scenarios$scenario_id, function(scenario_id) {
    do.call(rbind, lapply(PAIRWISE_N_VALUES, function(n) {
      data.frame(
        scenario_id = scenario_id,
        n = as.integer(n),
        replicate = seq_len(PAIRWISE_FULL_REPLICATES),
        stringsAsFactors = FALSE
      )
    }))
  }))
  full_map$stream_id <- seq_len(nrow(full_map))
  full_map$study_component <- "pairwise_full"
  full_map <- full_map[c("stream_id", "scenario_id", "n", "replicate", "study_component")]
  full_streams <- make_streams(nrow(full_map), MASTER_SEEDS[["pairwise"]])

  bootstrap_map <- do.call(rbind, lapply(c("S3", "S6"), function(scenario_id) {
    do.call(rbind, lapply(c(50L, 100L), function(n) {
      data.frame(
        scenario_id = scenario_id,
        n = n,
        outer_replicate = seq_len(1000L),
        bootstrap_draws = 999L,
        stringsAsFactors = FALSE
      )
    }))
  }))
  bootstrap_map$stream_id <- seq_len(nrow(bootstrap_map))
  bootstrap_map$study_component <- "focused_interval_bootstrap"
  bootstrap_map <- bootstrap_map[c(
    "stream_id", "scenario_id", "n", "outer_replicate", "bootstrap_draws", "study_component"
  )]
  bootstrap_streams <- make_streams(nrow(bootstrap_map), MASTER_SEEDS[["interval_bootstrap"]])

  rmse_map <- expand.grid(
    scenario_id = scenarios$scenario_id,
    n = PAIRWISE_N_VALUES,
    method = sort(c("proposed", "naive", "raw", "oracle_sample")),
    stringsAsFactors = FALSE
  )
  rmse_map <- rmse_map[order(rmse_map$scenario_id, rmse_map$n, rmse_map$method), , drop = FALSE]
  rownames(rmse_map) <- NULL
  rmse_map$stream_id <- seq_len(nrow(rmse_map))
  rmse_map$study_component <- "full_rmse_mcse"
  rmse_map <- rmse_map[c("stream_id", "scenario_id", "n", "method", "study_component")]
  rmse_streams <- make_streams(nrow(rmse_map), MASTER_SEEDS[["rmse_mcse"]])

  stream_dir <- file.path(root, "config", "rng_streams")
  dir.create(stream_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(
    pairwise_rds = file.path(stream_dir, "pairwise_full_streams.rds"),
    pairwise_csv = file.path(stream_dir, "pairwise_full_stream_map.csv"),
    bootstrap_rds = file.path(stream_dir, "interval_bootstrap_streams.rds"),
    bootstrap_csv = file.path(stream_dir, "interval_bootstrap_stream_map.csv"),
    rmse_rds = file.path(stream_dir, "rmse_full_streams.rds"),
    rmse_csv = file.path(stream_dir, "rmse_full_stream_map.csv")
  )
  atomic_save_rds(full_streams, files[["pairwise_rds"]])
  write.csv(full_map, files[["pairwise_csv"]], row.names = FALSE, quote = TRUE)
  atomic_save_rds(bootstrap_streams, files[["bootstrap_rds"]])
  write.csv(bootstrap_map, files[["bootstrap_csv"]], row.names = FALSE, quote = TRUE)
  atomic_save_rds(rmse_streams, files[["rmse_rds"]])
  write.csv(rmse_map, files[["rmse_csv"]], row.names = FALSE, quote = TRUE)

  stopifnot(
    length(readRDS(files[["pairwise_rds"]])) == 120000L,
    nrow(read.csv(files[["pairwise_csv"]])) == 120000L,
    length(readRDS(files[["bootstrap_rds"]])) == 4000L,
    nrow(read.csv(files[["bootstrap_csv"]])) == 4000L,
    length(readRDS(files[["rmse_rds"]])) == 96L,
    nrow(read.csv(files[["rmse_csv"]])) == 96L
  )
  data.frame(
    file = vapply(files, phase3_relative_path, character(1), root = root),
    sha256 = vapply(files, sha256_file, character(1)),
    bytes = vapply(files, function(path) file.info(path)$size, numeric(1)),
    master_seed = c(
      MASTER_SEEDS[["pairwise"]], NA, MASTER_SEEDS[["interval_bootstrap"]], NA,
      MASTER_SEEDS[["rmse_mcse"]], NA
    ),
    stringsAsFactors = FALSE
  )
}

canonical_bootstrap_tasks <- function(streams, mapping) {
  required <- c(
    "stream_id", "scenario_id", "n", "outer_replicate", "bootstrap_draws", "study_component"
  )
  if (!all(required %in% names(mapping))) stop("Invalid bootstrap RNG stream mapping schema", call. = FALSE)
  if (length(streams) != nrow(mapping)) stop("Bootstrap stream count does not match mapping", call. = FALSE)
  mapping <- mapping[order(mapping$scenario_id, mapping$n, mapping$outer_replicate), , drop = FALSE]
  rownames(mapping) <- NULL
  lapply(seq_len(nrow(mapping)), function(i) {
    list(
      stream_id = as.integer(mapping$stream_id[i]),
      scenario_id = mapping$scenario_id[i],
      n = as.integer(mapping$n[i]),
      outer_replicate = as.integer(mapping$outer_replicate[i]),
      bootstrap_draws = as.integer(mapping$bootstrap_draws[i]),
      seed = streams[[mapping$stream_id[i]]]
    )
  })
}

.metric_mcse_bootstrap <- function(estimate, truth, stream, resamples = 2000L) {
  estimate <- estimate[is.finite(estimate)]
  if (length(estimate) < 2L) {
    return(c(empirical_sd = NA_real_, rmse = NA_real_, median_absolute_error = NA_real_))
  }
  assign(".Random.seed", stream, envir = .GlobalEnv)
  B <- length(estimate)
  values <- matrix(NA_real_, nrow = resamples, ncol = 3L)
  for (i in seq_len(resamples)) {
    draw <- estimate[sample.int(B, B, replace = TRUE)]
    error <- draw - truth
    values[i, ] <- c(stats::sd(draw), sqrt(mean(error^2)), stats::median(abs(error)))
  }
  stats::setNames(apply(values, 2L, stats::sd), c(
    "empirical_sd", "rmse", "median_absolute_error"
  ))
}

summarise_pairwise_full <- function(results, rmse_streams, rmse_mapping) {
  validate_pairwise_pilot_output(results, 120000L)
  methods <- c(
    proposed = "rho_proposed", naive = "rho_naive",
    raw = "rho_raw", oracle_sample = "rho_oracle_sample"
  )
  groups <- unique(results[c("scenario_id", "n")])
  groups <- groups[order(groups$scenario_id, groups$n), , drop = FALSE]
  rows <- vector("list", nrow(groups) * length(methods))
  z <- 0L
  for (g in seq_len(nrow(groups))) {
    cell <- results[
      results$scenario_id == groups$scenario_id[g] & results$n == groups$n[g], , drop = FALSE
    ]
    truth <- unique(cell$rho_true)
    if (length(truth) != 1L) stop("Non-unique analytic truth within full pairwise cell", call. = FALSE)
    for (method in names(methods)) {
      z <- z + 1L
      estimate <- cell[[methods[[method]]]]
      valid <- is.finite(estimate)
      error <- estimate[valid] - truth
      stream_row <- rmse_mapping$scenario_id == groups$scenario_id[g] &
        rmse_mapping$n == groups$n[g] & rmse_mapping$method == method
      if (sum(stream_row) != 1L) stop("Missing full-study metric-MCSE stream mapping", call. = FALSE)
      mcse <- .metric_mcse_bootstrap(
        estimate[valid], truth, rmse_streams[[rmse_mapping$stream_id[stream_row]]]
      )
      proposed <- identical(method, "proposed")
      coverage <- if (proposed) mean(cell$if_cover) else NA_real_
      lower_miss <- if (proposed) mean(cell$if_lower > truth) else NA_real_
      upper_miss <- if (proposed) mean(cell$if_upper < truth) else NA_real_
      fisher_sd <- if (proposed) stats::sd(atanh(estimate[valid])) else NA_real_
      mean_fisher_se <- if (proposed) mean(cell$if_se_fisher) else NA_real_
      rows[[z]] <- data.frame(
        scenario_id = groups$scenario_id[g], n = groups$n[g], method = method,
        B_total = nrow(cell), B_valid = sum(valid), truth = truth,
        mean_estimate = mean(estimate[valid]), bias = mean(error),
        mcse_bias = stats::sd(error) / sqrt(length(error)),
        empirical_sd = stats::sd(estimate[valid]),
        mcse_empirical_sd = unname(mcse[["empirical_sd"]]),
        rmse = sqrt(mean(error^2)), mcse_rmse = unname(mcse[["rmse"]]),
        median_absolute_error = stats::median(abs(error)),
        mcse_median_absolute_error = unname(mcse[["median_absolute_error"]]),
        coverage = coverage,
        mcse_coverage = if (proposed) sqrt(coverage * (1 - coverage) / nrow(cell)) else NA_real_,
        mean_interval_length = if (proposed) mean(cell$if_length) else NA_real_,
        mcse_mean_interval_length = if (proposed) stats::sd(cell$if_length) / sqrt(nrow(cell)) else NA_real_,
        lower_noncoverage = lower_miss,
        mcse_lower_noncoverage = if (proposed) sqrt(lower_miss * (1 - lower_miss) / nrow(cell)) else NA_real_,
        upper_noncoverage = upper_miss,
        mcse_upper_noncoverage = if (proposed) sqrt(upper_miss * (1 - upper_miss) / nrow(cell)) else NA_real_,
        mean_estimated_se_fisher = mean_fisher_se,
        mcse_mean_estimated_se_fisher = if (proposed) stats::sd(cell$if_se_fisher) / sqrt(nrow(cell)) else NA_real_,
        empirical_sd_fisher = fisher_sd,
        se_fisher_to_empirical_fisher_sd_ratio = if (proposed) mean_fisher_se / fisher_sd else NA_real_,
        projection_rate = mean(cell$projection_active),
        negative_variance_rate = mean(cell$negative_raw_variance),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarise_pairwise_diagnostics <- function(results) {
  validate_pairwise_pilot_output(results, 120000L)
  groups <- unique(results[c("scenario_id", "n")])
  groups <- groups[order(groups$scenario_id, groups$n), , drop = FALSE]
  do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    cell <- results[
      results$scenario_id == groups$scenario_id[i] & results$n == groups$n[i], , drop = FALSE
    ]
    data.frame(
      scenario_id = groups$scenario_id[i], n = groups$n[i], replicates = nrow(cell),
      failures = sum(!cell$success), warnings = sum(cell$warning_count),
      unexplained_nonfinite = sum(cell$nonfinite_count),
      raw_negative_variance_count = sum(cell$negative_raw_variance),
      raw_negative_variance_rate = mean(cell$negative_raw_variance),
      projection_count = sum(cell$projection_active), projection_rate = mean(cell$projection_active),
      floor_count = sum(cell$floor_active), floor_rate = mean(cell$floor_active),
      min_eigen_min = min(cell$min_eigen_raw),
      min_eigen_q05 = unname(stats::quantile(cell$min_eigen_raw, 0.05, type = 7)),
      min_eigen_median = stats::median(cell$min_eigen_raw),
      min_eigen_q95 = unname(stats::quantile(cell$min_eigen_raw, 0.95, type = 7)),
      mean_M1 = mean(cell$mean_M1), mean_M2 = mean(cell$mean_M2),
      median_M1 = stats::median(cell$median_M1), median_M2 = stats::median(cell$median_M2),
      mean_cor_M1_M2 = mean(cell$cor_M1_M2),
      mean_cor_M1_P1 = mean(cell$cor_M1_P1), mean_cor_M2_P2 = mean(cell$cor_M2_P2),
      stringsAsFactors = FALSE
    )
  }))
}

.bootstrap_result_template <- function(task) {
  data.frame(
    scenario_id = as.character(task$scenario_id), n = as.integer(task$n),
    outer_replicate = as.integer(task$outer_replicate), stream_id = as.integer(task$stream_id),
    outer_seed_hash = "", bootstrap_seed_hash = "", rho_true = NA_real_,
    rho_proposed = NA_real_, analytic_lower = NA_real_, analytic_upper = NA_real_,
    analytic_cover = NA, analytic_length = NA_real_, analytic_lower_miss = NA,
    analytic_upper_miss = NA, bootstrap_lower = NA_real_, bootstrap_upper = NA_real_,
    bootstrap_cover = NA, bootstrap_length = NA_real_, bootstrap_lower_miss = NA,
    bootstrap_upper_miss = NA, bootstrap_draws_target = as.integer(task$bootstrap_draws),
    bootstrap_draws_attempted = 0L, bootstrap_draws_valid = 0L,
    bootstrap_failure_count = 0L, bootstrap_failure_reasons = "",
    bootstrap_ci_valid = FALSE, quantile_type = 8L,
    projection_active = NA, floor_active = NA, near_boundary = NA,
    bootstrap_projection_count = 0L, bootstrap_floor_count = 0L,
    bootstrap_boundary_count = 0L, warning_count = 0L, warning_messages = "",
    bootstrap_warning_count = 0L, success = FALSE, failure_message = "",
    probability_bounds_valid = FALSE, denominator_valid = FALSE, count_valid = FALSE,
    denominator_mechanism_valid = FALSE, unweighted_valid = FALSE,
    rho_truth_analytic = FALSE, runtime_seconds = NA_real_, stringsAsFactors = FALSE
  )
}

.collapse_reason_counts <- function(reasons) {
  if (!length(reasons)) return("")
  counts <- sort(table(reasons), decreasing = TRUE)
  paste(paste(names(counts), as.integer(counts), sep = "="), collapse = " | ")
}

focused_bootstrap_replicate <- function(task) {
  out <- .bootstrap_result_template(task)
  warnings <- character()
  bootstrap_warnings <- character()
  failure_reasons <- character()
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    withCallingHandlers({
      assign(".Random.seed", task$seed, envir = .GlobalEnv)
      out$outer_seed_hash <- digest::digest(task$seed, algo = "sha256", serialize = TRUE)
      dat <- generate_pairwise_data(task$n, task$scenario_id)
      out$bootstrap_seed_hash <- digest::digest(.Random.seed, algo = "sha256", serialize = TRUE)
      fit <- fit_latent_binomial_cov(dat$X, dat$M, keep_influence = TRUE)
      ci <- latent_correlation_ci(fit, 1L, 2L, level = 0.95, scale = "fisher")
      truth <- dat$truth$rho

      component <- dat$denominator_components
      mechanism_ok <- if (dat$scenario$denominator_id == "D1") {
        identical(as.numeric(dat$M[, 1L]), as.numeric(2L + component$N1)) &&
          identical(as.numeric(dat$M[, 2L]), as.numeric(2L + component$N2))
      } else if (dat$scenario$denominator_id == "D2") {
        identical(as.numeric(dat$M[, 1L]), as.numeric(2L + component$N1)) &&
          identical(as.numeric(dat$M[, 2L]), as.numeric(2L + component$N2))
      } else {
        identical(as.numeric(dat$M[, 1L]), as.numeric(2L + component$C + component$E1)) &&
          identical(as.numeric(dat$M[, 2L]), as.numeric(2L + component$C + component$E2)) &&
          isTRUE(component$shared_C_used_once)
      }
      out$probability_bounds_valid <- all(is.finite(dat$P)) && all(dat$P >= 0 & dat$P <= 1)
      out$denominator_valid <- all(dat$M >= 2) && all(abs(dat$M - round(dat$M)) <= INTEGER_TOL)
      out$count_valid <- all(abs(dat$X - round(dat$X)) <= INTEGER_TOL) &&
        all(dat$X >= 0 & dat$X <= dat$M)
      out$denominator_mechanism_valid <- mechanism_ok
      out$unweighted_valid <- isTRUE(fit$validation_summary$unweighted_rows)
      out$rho_truth_analytic <- identical(unname(truth), unname(pairwise_truth(task$scenario_id)$rho))
      out$rho_true <- truth
      out$rho_proposed <- fit$R_latent[1L, 2L]
      out$analytic_lower <- ci$lower
      out$analytic_upper <- ci$upper
      out$analytic_cover <- ci$lower <= truth && ci$upper >= truth
      out$analytic_length <- ci$upper - ci$lower
      out$analytic_lower_miss <- ci$lower > truth
      out$analytic_upper_miss <- ci$upper < truth
      out$projection_active <- fit$projection_active
      out$floor_active <- fit$floor_active
      out$near_boundary <- abs(out$rho_proposed) > 0.95

      z_boot <- rep(NA_real_, task$bootstrap_draws)
      for (b in seq_len(task$bootstrap_draws)) {
        out$bootstrap_draws_attempted <- b
        index <- sample.int(task$n, task$n, replace = TRUE)
        boot_fit <- tryCatch(
          withCallingHandlers(
            fit_latent_binomial_cov(dat$X[index, , drop = FALSE], dat$M[index, , drop = FALSE],
                                    keep_influence = FALSE),
            warning = function(w) {
              bootstrap_warnings <<- c(bootstrap_warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          ),
          error = function(e) e
        )
        if (inherits(boot_fit, "error")) {
          failure_reasons <- c(failure_reasons, paste0("fit_error: ", conditionMessage(boot_fit)))
          next
        }
        rho_boot <- boot_fit$R_latent[1L, 2L]
        out$bootstrap_projection_count <- out$bootstrap_projection_count + as.integer(boot_fit$projection_active)
        out$bootstrap_floor_count <- out$bootstrap_floor_count + as.integer(boot_fit$floor_active)
        out$bootstrap_boundary_count <- out$bootstrap_boundary_count + as.integer(abs(rho_boot) > 0.95)
        if (!is.finite(rho_boot) || abs(rho_boot) >= 1) {
          failure_reasons <- c(failure_reasons, "nonfinite_or_boundary_fisher_transform")
          next
        }
        z_boot[b] <- atanh(rho_boot)
      }
      valid <- is.finite(z_boot)
      out$bootstrap_draws_valid <- sum(valid)
      out$bootstrap_failure_count <- sum(!valid)
      out$bootstrap_failure_reasons <- .collapse_reason_counts(failure_reasons)
      out$bootstrap_warning_count <- length(bootstrap_warnings)
      if (all(valid) && length(z_boot) == task$bootstrap_draws) {
        z_limits <- stats::quantile(z_boot, probs = c(0.025, 0.975), type = 8, names = FALSE)
        limits <- tanh(z_limits)
        out$bootstrap_lower <- limits[1L]
        out$bootstrap_upper <- limits[2L]
        out$bootstrap_cover <- limits[1L] <= truth && limits[2L] >= truth
        out$bootstrap_length <- limits[2L] - limits[1L]
        out$bootstrap_lower_miss <- limits[1L] > truth
        out$bootstrap_upper_miss <- limits[2L] < truth
        out$bootstrap_ci_valid <- all(is.finite(limits)) && limits[1L] <= limits[2L]
      }
      out$success <- all(c(
        out$probability_bounds_valid, out$denominator_valid, out$count_valid,
        out$denominator_mechanism_valid, out$unweighted_valid, out$rho_truth_analytic,
        is.finite(out$rho_proposed), is.finite(out$analytic_lower), is.finite(out$analytic_upper),
        out$bootstrap_draws_attempted == task$bootstrap_draws,
        out$bootstrap_draws_valid == task$bootstrap_draws, out$bootstrap_ci_valid
      ))
      out
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      out$failure_message <- conditionMessage(e)
      out
    }
  )
  result$warning_count <- length(warnings)
  result$warning_messages <- paste(unique(warnings), collapse = " | ")
  result$runtime_seconds <- proc.time()[["elapsed"]] - started
  result
}

canonical_bootstrap_results <- function(x) {
  if (!nrow(x)) return(x)
  x <- x[order(x$scenario_id, x$n, x$outer_replicate), , drop = FALSE]
  rownames(x) <- NULL
  x
}

execute_bootstrap_tasks <- function(tasks, cluster = NULL, load_balanced = TRUE) {
  if (!length(tasks)) return(data.frame())
  rows <- if (is.null(cluster)) {
    lapply(tasks, focused_bootstrap_replicate)
  } else if (isTRUE(load_balanced)) {
    parallel::parLapplyLB(cluster, tasks, function(task) focused_bootstrap_replicate(task))
  } else {
    parallel::parLapply(cluster, tasks, function(task) focused_bootstrap_replicate(task))
  }
  canonical_bootstrap_results(do.call(rbind, rows))
}

.bootstrap_task_key <- function(task) paste(task$scenario_id, task$n, task$outer_replicate, sep = "|")
.bootstrap_result_keys <- function(results) {
  if (!nrow(results)) character() else paste(results$scenario_id, results$n, results$outer_replicate, sep = "|")
}

run_checkpointed_bootstrap_tasks <- function(tasks, checkpoint_path, final_path, cluster = NULL,
                                             checkpoint_every = 10L, load_balanced = TRUE) {
  existing <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else data.frame()
  completed <- .bootstrap_result_keys(existing)
  pending <- tasks[!vapply(tasks, .bootstrap_task_key, character(1)) %in% completed]
  if (length(pending)) {
    starts <- seq.int(1L, length(pending), by = as.integer(checkpoint_every))
    for (start in starts) {
      chunk <- pending[start:min(start + checkpoint_every - 1L, length(pending))]
      new <- execute_bootstrap_tasks(chunk, cluster, load_balanced)
      existing <- if (nrow(existing)) rbind(existing, new) else new
      existing <- canonical_bootstrap_results(existing)
      atomic_save_rds(existing, checkpoint_path)
    }
  }
  all_keys <- vapply(tasks, .bootstrap_task_key, character(1))
  complete <- setequal(.bootstrap_result_keys(existing), all_keys) && nrow(existing) == length(tasks)
  if (complete) atomic_save_rds(existing, final_path)
  list(
    results = existing, complete = complete, completed = nrow(existing), planned = length(tasks),
    checkpoint_path = checkpoint_path, final_path = final_path
  )
}

bootstrap_required_columns <- function() c(
  "scenario_id", "n", "outer_replicate", "stream_id", "outer_seed_hash",
  "bootstrap_seed_hash", "rho_true", "rho_proposed", "analytic_lower", "analytic_upper",
  "analytic_cover", "analytic_length", "analytic_lower_miss", "analytic_upper_miss",
  "bootstrap_lower", "bootstrap_upper", "bootstrap_cover", "bootstrap_length",
  "bootstrap_lower_miss", "bootstrap_upper_miss", "bootstrap_draws_target",
  "bootstrap_draws_attempted", "bootstrap_draws_valid", "bootstrap_failure_count",
  "bootstrap_failure_reasons", "bootstrap_ci_valid", "quantile_type", "projection_active",
  "floor_active", "near_boundary", "bootstrap_projection_count", "bootstrap_floor_count",
  "bootstrap_boundary_count", "warning_count", "warning_messages", "bootstrap_warning_count",
  "success", "failure_message", "probability_bounds_valid", "denominator_valid", "count_valid",
  "denominator_mechanism_valid", "unweighted_valid", "rho_truth_analytic", "runtime_seconds"
)

validate_bootstrap_output <- function(x, expected_rows = 4000L) {
  missing <- setdiff(bootstrap_required_columns(), names(x))
  if (length(missing)) stop("Bootstrap output missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(x) != expected_rows) stop("Unexpected bootstrap outer-dataset count", call. = FALSE)
  keys <- .bootstrap_result_keys(x)
  if (anyDuplicated(keys)) stop("Duplicate bootstrap outer task keys", call. = FALSE)
  if (any(!x$success)) stop("Focused bootstrap contains failed outer tasks", call. = FALSE)
  if (any(x$bootstrap_draws_target != 999L) || any(x$bootstrap_draws_attempted != 999L)) {
    stop("Focused bootstrap did not attempt 999 draws for every outer dataset", call. = FALSE)
  }
  if (any(x$bootstrap_draws_valid != 999L) || any(x$bootstrap_failure_count != 0L)) {
    stop("Focused bootstrap contains invalid or failed inner draws", call. = FALSE)
  }
  if (any(x$quantile_type != 8L) || any(!x$bootstrap_ci_valid)) {
    stop("Focused bootstrap interval construction failed its schema", call. = FALSE)
  }
  contract <- c(
    "probability_bounds_valid", "denominator_valid", "count_valid",
    "denominator_mechanism_valid", "unweighted_valid", "rho_truth_analytic"
  )
  if (!all(vapply(x[contract], all, logical(1)))) stop("Bootstrap generator contract failed", call. = FALSE)
  invisible(TRUE)
}

summarise_focused_bootstrap <- function(results) {
  validate_bootstrap_output(results)
  groups <- unique(results[c("scenario_id", "n")])
  groups <- groups[order(groups$scenario_id, groups$n), , drop = FALSE]
  do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    cell <- results[
      results$scenario_id == groups$scenario_id[i] & results$n == groups$n[i], , drop = FALSE
    ]
    analytic_coverage <- mean(cell$analytic_cover)
    bootstrap_coverage <- mean(cell$bootstrap_cover)
    data.frame(
      scenario_id = groups$scenario_id[i], n = groups$n[i], outer_datasets = nrow(cell),
      all_999_valid = sum(cell$bootstrap_draws_valid == 999L),
      total_bootstrap_draws_attempted = sum(cell$bootstrap_draws_attempted),
      total_bootstrap_draws_valid = sum(cell$bootstrap_draws_valid),
      total_bootstrap_failures = sum(cell$bootstrap_failure_count),
      minimum_valid_bootstrap_draws = min(cell$bootstrap_draws_valid),
      median_valid_bootstrap_draws = stats::median(cell$bootstrap_draws_valid),
      maximum_valid_bootstrap_draws = max(cell$bootstrap_draws_valid),
      analytic_coverage = analytic_coverage,
      analytic_mcse_coverage = sqrt(analytic_coverage * (1 - analytic_coverage) / nrow(cell)),
      bootstrap_coverage = bootstrap_coverage,
      bootstrap_mcse_coverage = sqrt(bootstrap_coverage * (1 - bootstrap_coverage) / nrow(cell)),
      coverage_difference_bootstrap_minus_analytic = bootstrap_coverage - analytic_coverage,
      analytic_mean_length = mean(cell$analytic_length),
      analytic_mcse_mean_length = stats::sd(cell$analytic_length) / sqrt(nrow(cell)),
      bootstrap_mean_length = mean(cell$bootstrap_length),
      bootstrap_mcse_mean_length = stats::sd(cell$bootstrap_length) / sqrt(nrow(cell)),
      analytic_lower_noncoverage = mean(cell$analytic_lower_miss),
      analytic_upper_noncoverage = mean(cell$analytic_upper_miss),
      bootstrap_lower_noncoverage = mean(cell$bootstrap_lower_miss),
      bootstrap_upper_noncoverage = mean(cell$bootstrap_upper_miss),
      outer_projection_rate = mean(cell$projection_active),
      outer_floor_rate = mean(cell$floor_active),
      outer_boundary_rate = mean(cell$near_boundary),
      bootstrap_projection_rate = sum(cell$bootstrap_projection_count) / sum(cell$bootstrap_draws_attempted),
      bootstrap_floor_rate = sum(cell$bootstrap_floor_count) / sum(cell$bootstrap_draws_attempted),
      bootstrap_boundary_rate = sum(cell$bootstrap_boundary_count) / sum(cell$bootstrap_draws_attempted),
      analytic_warning_count = sum(cell$warning_count),
      bootstrap_warning_count = sum(cell$bootstrap_warning_count),
      ci_construction_failures = sum(!cell$bootstrap_ci_valid),
      outer_runtime_seconds = sum(cell$runtime_seconds),
      stringsAsFactors = FALSE
    )
  }))
}
