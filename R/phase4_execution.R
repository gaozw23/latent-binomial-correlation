.phase4_project_root <- function() {
  normalizePath(getOption("latent_binomial_project_root", "."), winslash = "/", mustWork = TRUE)
}

generate_phase4_stream_files <- function(root = .phase4_project_root()) {
  mapping <- do.call(rbind, lapply(MATRIX_N_VALUES, function(n) {
    data.frame(
      n = as.integer(n), replicate = seq_len(MATRIX_FULL_REPLICATES),
      stringsAsFactors = FALSE
    )
  }))
  mapping$stream_id <- seq_len(nrow(mapping))
  mapping$study_component <- "matrix_full"
  mapping <- mapping[c("stream_id", "n", "replicate", "study_component")]
  streams <- make_streams(nrow(mapping), MASTER_SEEDS[["matrix"]])

  stream_dir <- file.path(root, "config", "rng_streams")
  dir.create(stream_dir, recursive = TRUE, showWarnings = FALSE)
  rds_path <- file.path(stream_dir, "matrix_full_streams.rds")
  map_path <- file.path(stream_dir, "matrix_full_stream_map.csv")
  atomic_save_rds(streams, rds_path)
  write.csv(mapping, map_path, row.names = FALSE, quote = TRUE)
  stopifnot(length(readRDS(rds_path)) == 15000L, nrow(read.csv(map_path)) == 15000L)
  data.frame(
    file = c("config/rng_streams/matrix_full_streams.rds",
             "config/rng_streams/matrix_full_stream_map.csv"),
    sha256 = c(sha256_file(rds_path), sha256_file(map_path)),
    bytes = c(file.info(rds_path)$size, file.info(map_path)$size),
    master_seed = c(MASTER_SEEDS[["matrix"]], NA), stringsAsFactors = FALSE
  )
}

canonical_matrix_full_tasks <- function(streams, mapping) {
  canonical_matrix_tasks(streams, mapping)
}

.matrix_list_column <- function(value) I(list(value))

.matrix_full_result_template <- function(task) {
  empty <- matrix(NA_real_, 5L, 5L)
  data.frame(
    n = as.integer(task$n), replicate = as.integer(task$replicate),
    stream_id = as.integer(task$stream_id),
    Sigma_raw = .matrix_list_column(empty), Sigma_psd = .matrix_list_column(empty),
    Sigma_stable = .matrix_list_column(empty), R_estimate = .matrix_list_column(empty),
    Partial_estimate = .matrix_list_column(empty), Sigma_true = .matrix_list_column(empty),
    R_true = .matrix_list_column(empty), Theta_true = .matrix_list_column(empty),
    Partial_true = .matrix_list_column(empty),
    loss_raw_fro = NA_real_, loss_psd_fro = NA_real_, loss_stable_fro = NA_real_,
    loss_raw_operator = NA_real_, loss_psd_operator = NA_real_, loss_stable_operator = NA_real_,
    loss_ratio_psd_raw = NA_real_, projection_improves_loss = NA,
    max_abs_R_error = NA_real_, fro_R_error = NA_real_, rmse_offdiag_R = NA_real_,
    max_abs_partial_error = NA_real_, fro_partial_error = NA_real_,
    rmse_offdiag_partial = NA_real_, min_eigen_raw = NA_real_,
    projection_active = NA, floor_active = NA, projection_inequality_pass = NA,
    projection_inequality_margin = NA_real_, success = FALSE,
    warning_count = 0L, warning_messages = "", nonfinite_count = NA_integer_,
    failure_message = "", probability_bounds_valid = FALSE,
    denominator_valid = FALSE, count_valid = FALSE, denominator_mechanism_valid = FALSE,
    W_valid = FALSE, truth_valid = FALSE, unweighted_valid = FALSE,
    stringsAsFactors = FALSE
  )
}

matrix_full_replicate <- function(task) {
  out <- .matrix_full_result_template(task)
  warnings <- character()
  result <- tryCatch(
    withCallingHandlers({
      assign(".Random.seed", task$seed, envir = .GlobalEnv)
      dat <- generate_matrix_data(task$n)
      fit <- fit_latent_binomial_cov(dat$X, dat$M, keep_influence = FALSE)
      truth <- dat$truth

      sigma_expected <- truth$sigma2_B * truth$W %*% t(truth$W)
      W_ok <- identical(dim(truth$W), c(5L, 8L)) && all(truth$W >= 0) &&
        all(rowSums(truth$W) == 1)
      truth_ok <- isTRUE(all.equal(truth$Sigma, sigma_expected, tolerance = 0)) &&
        min(eigen(truth$Sigma, symmetric = TRUE, only.values = TRUE)$values) > 0 &&
        isTRUE(all.equal(truth$R, cov2cor(truth$Sigma), tolerance = 1e-15)) &&
        isTRUE(all.equal(truth$Theta, solve(truth$Sigma), tolerance = 1e-12))
      probability_ok <- all(is.finite(dat$P)) && all(dat$P >= 0 & dat$P <= 1)
      denominator_ok <- all(dat$M >= 2) && all(abs(dat$M - round(dat$M)) <= INTEGER_TOL)
      count_ok <- all(abs(dat$X - round(dat$X)) <= INTEGER_TOL) &&
        all(dat$X >= 0 & dat$X <= dat$M)
      mechanism_ok <- identical(unname(dat$M), unname(2L + dat$C + dat$EM))

      raw_error <- symmetrise_matrix(fit$Sigma_raw - truth$Sigma)
      psd_error <- symmetrise_matrix(fit$Sigma_psd - truth$Sigma)
      stable_error <- symmetrise_matrix(fit$Sigma_stable - truth$Sigma)
      raw_fro <- norm(raw_error, "F")
      psd_fro <- norm(psd_error, "F")
      stable_fro <- norm(stable_error, "F")
      theorem <- frobenius_projection_check(fit$Sigma_raw, fit$Sigma_psd, truth$Sigma)
      R_error <- fit$R_latent - truth$R
      partial_error <- fit$Partial - truth$Partial
      numeric_required <- c(
        fit$Sigma_raw, fit$Sigma_psd, fit$Sigma_stable, fit$R_latent, fit$Partial,
        truth$Sigma, truth$R, truth$Theta, truth$Partial,
        raw_fro, psd_fro, stable_fro, theorem$margin, R_error, partial_error
      )
      nonfinite <- sum(!is.finite(numeric_required))

      out$Sigma_raw[[1L]] <- fit$Sigma_raw
      out$Sigma_psd[[1L]] <- fit$Sigma_psd
      out$Sigma_stable[[1L]] <- fit$Sigma_stable
      out$R_estimate[[1L]] <- fit$R_latent
      out$Partial_estimate[[1L]] <- fit$Partial
      out$Sigma_true[[1L]] <- truth$Sigma
      out$R_true[[1L]] <- truth$R
      out$Theta_true[[1L]] <- truth$Theta
      out$Partial_true[[1L]] <- truth$Partial
      out$loss_raw_fro <- raw_fro
      out$loss_psd_fro <- psd_fro
      out$loss_stable_fro <- stable_fro
      out$loss_raw_operator <- max(abs(eigen(raw_error, symmetric = TRUE, only.values = TRUE)$values))
      out$loss_psd_operator <- max(abs(eigen(psd_error, symmetric = TRUE, only.values = TRUE)$values))
      out$loss_stable_operator <- max(abs(eigen(stable_error, symmetric = TRUE, only.values = TRUE)$values))
      out$loss_ratio_psd_raw <- if (raw_fro == 0) 0 else psd_fro / raw_fro
      out$projection_improves_loss <-
        (raw_fro - psd_fro) > PROJECTION_LOSS_TOL * (1 + raw_fro)
      out$max_abs_R_error <- max(abs(R_error))
      out$fro_R_error <- norm(R_error, "F")
      out$rmse_offdiag_R <- sqrt(mean(R_error[upper.tri(R_error)]^2))
      out$max_abs_partial_error <- max(abs(partial_error))
      out$fro_partial_error <- norm(partial_error, "F")
      out$rmse_offdiag_partial <- sqrt(mean(partial_error[upper.tri(partial_error)]^2))
      out$min_eigen_raw <- fit$min_eigen_raw
      out$projection_active <- fit$projection_active
      out$floor_active <- fit$floor_active
      out$projection_inequality_pass <- theorem$pass
      out$projection_inequality_margin <- theorem$margin
      out$nonfinite_count <- as.integer(nonfinite)
      out$probability_bounds_valid <- probability_ok
      out$denominator_valid <- denominator_ok
      out$count_valid <- count_ok
      out$denominator_mechanism_valid <- mechanism_ok
      out$W_valid <- W_ok
      out$truth_valid <- truth_ok
      out$unweighted_valid <- isTRUE(fit$validation_summary$unweighted_rows)
      out$success <- all(c(
        probability_ok, denominator_ok, count_ok, mechanism_ok, W_ok, truth_ok,
        out$unweighted_valid, theorem$pass, nonfinite == 0L
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
  result
}

canonical_matrix_full_results <- function(x) {
  if (!nrow(x)) return(x)
  x <- x[order(x$n, x$replicate), , drop = FALSE]
  rownames(x) <- NULL
  x
}

execute_matrix_full_tasks <- function(tasks, cluster = NULL, load_balanced = TRUE) {
  if (!length(tasks)) return(data.frame())
  rows <- if (is.null(cluster)) {
    lapply(tasks, matrix_full_replicate)
  } else if (isTRUE(load_balanced)) {
    parallel::parLapplyLB(cluster, tasks, function(task) matrix_full_replicate(task))
  } else {
    parallel::parLapply(cluster, tasks, function(task) matrix_full_replicate(task))
  }
  canonical_matrix_full_results(do.call(rbind, rows))
}

.matrix_full_task_key <- function(task) paste(task$n, task$replicate, sep = "|")
.matrix_full_result_keys <- function(results) {
  if (!nrow(results)) character() else paste(results$n, results$replicate, sep = "|")
}

run_checkpointed_matrix_full_tasks <- function(tasks, checkpoint_path, final_path,
                                               cluster = NULL, checkpoint_every = 500L,
                                               load_balanced = TRUE) {
  existing <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else data.frame()
  completed <- .matrix_full_result_keys(existing)
  pending <- tasks[!vapply(tasks, .matrix_full_task_key, character(1)) %in% completed]
  if (length(pending)) {
    starts <- seq.int(1L, length(pending), by = as.integer(checkpoint_every))
    for (start in starts) {
      chunk <- pending[start:min(start + checkpoint_every - 1L, length(pending))]
      new <- execute_matrix_full_tasks(chunk, cluster, load_balanced)
      existing <- if (nrow(existing)) rbind(existing, new) else new
      existing <- canonical_matrix_full_results(existing)
      atomic_save_rds(existing, checkpoint_path)
    }
  }
  all_keys <- vapply(tasks, .matrix_full_task_key, character(1))
  complete <- setequal(.matrix_full_result_keys(existing), all_keys) && nrow(existing) == length(tasks)
  if (complete) atomic_save_rds(existing, final_path)
  list(results = existing, complete = complete, completed = nrow(existing), planned = length(tasks))
}

matrix_full_required_columns <- function() c(
  "n", "replicate", "stream_id", "Sigma_raw", "Sigma_psd", "Sigma_stable",
  "R_estimate", "Partial_estimate", "Sigma_true", "R_true", "Theta_true", "Partial_true",
  "loss_raw_fro", "loss_psd_fro", "loss_stable_fro", "loss_raw_operator",
  "loss_psd_operator", "loss_stable_operator", "loss_ratio_psd_raw",
  "projection_improves_loss", "max_abs_R_error", "fro_R_error", "rmse_offdiag_R",
  "max_abs_partial_error", "fro_partial_error", "rmse_offdiag_partial", "min_eigen_raw",
  "projection_active", "floor_active", "projection_inequality_pass",
  "projection_inequality_margin", "success", "warning_count", "warning_messages",
  "nonfinite_count", "failure_message", "probability_bounds_valid", "denominator_valid",
  "count_valid", "denominator_mechanism_valid", "W_valid", "truth_valid", "unweighted_valid"
)

validate_matrix_full_output <- function(x, expected_rows = 15000L) {
  missing <- setdiff(matrix_full_required_columns(), names(x))
  if (length(missing)) stop("Full matrix output missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(x) != expected_rows) stop("Unexpected full matrix row count", call. = FALSE)
  if (anyDuplicated(.matrix_full_result_keys(x))) stop("Duplicate full matrix task key", call. = FALSE)
  if (any(!x$success)) stop("Full matrix output contains failed tasks", call. = FALSE)
  if (sum(x$nonfinite_count) != 0L) stop("Full matrix output contains non-finite required results", call. = FALSE)
  if (any(!x$projection_inequality_pass)) stop("Full matrix PSD theorem violation", call. = FALSE)
  contract <- c(
    "probability_bounds_valid", "denominator_valid", "count_valid",
    "denominator_mechanism_valid", "W_valid", "truth_valid", "unweighted_valid"
  )
  if (!all(vapply(x[contract], all, logical(1)))) stop("Full matrix generator/truth contract failed", call. = FALSE)
  matrix_columns <- c(
    "Sigma_raw", "Sigma_psd", "Sigma_stable", "R_estimate", "Partial_estimate",
    "Sigma_true", "R_true", "Theta_true", "Partial_true"
  )
  for (column in matrix_columns) {
    if (length(x[[column]]) != nrow(x) ||
        !all(vapply(x[[column]], function(A) identical(dim(A), c(5L, 5L)) && all(is.finite(A)), logical(1)))) {
      stop("Invalid matrix list column: ", column, call. = FALSE)
    }
  }
  invisible(TRUE)
}

.mean_mcse <- function(x) c(mean = mean(x), mcse = sd(x) / sqrt(length(x)))
.rms_mcse <- function(x) {
  rms <- sqrt(mean(x^2))
  c(rms = rms, mcse = if (rms == 0) 0 else sd(x^2) / (2 * rms * sqrt(length(x))))
}

summarise_matrix_full <- function(results) {
  validate_matrix_full_output(results)
  do.call(rbind, lapply(sort(unique(results$n)), function(n) {
    x <- results[results$n == n, , drop = FALSE]
    raw_mean <- .mean_mcse(x$loss_raw_fro)
    psd_mean <- .mean_mcse(x$loss_psd_fro)
    stable_mean <- .mean_mcse(x$loss_stable_fro)
    raw_rms <- .rms_mcse(x$loss_raw_fro)
    psd_rms <- .rms_mcse(x$loss_psd_fro)
    stable_rms <- .rms_mcse(x$loss_stable_fro)
    ratio <- .mean_mcse(x$loss_ratio_psd_raw)
    R_max <- .mean_mcse(x$max_abs_R_error)
    R_fro <- .mean_mcse(x$fro_R_error)
    R_rmse <- .mean_mcse(x$rmse_offdiag_R)
    P_max <- .mean_mcse(x$max_abs_partial_error)
    P_fro <- .mean_mcse(x$fro_partial_error)
    P_rmse <- .mean_mcse(x$rmse_offdiag_partial)
    active <- x$projection_active | x$floor_active
    data.frame(
      n = n, replicates = nrow(x), failures = sum(!x$success),
      warning_count = sum(x$warning_count), unexplained_nonfinite = sum(x$nonfinite_count),
      projection_rate = mean(x$projection_active), floor_rate = mean(x$floor_active),
      min_eigen_min = min(x$min_eigen_raw),
      min_eigen_q05 = unname(quantile(x$min_eigen_raw, 0.05, type = 7)),
      min_eigen_median = median(x$min_eigen_raw),
      min_eigen_q95 = unname(quantile(x$min_eigen_raw, 0.95, type = 7)),
      raw_fro_mean = raw_mean[["mean"]], raw_fro_mean_mcse = raw_mean[["mcse"]],
      raw_fro_rms = raw_rms[["rms"]], raw_fro_rms_mcse = raw_rms[["mcse"]],
      raw_fro_q05 = unname(quantile(x$loss_raw_fro, 0.05, type = 7)),
      raw_fro_median = median(x$loss_raw_fro),
      raw_fro_q95 = unname(quantile(x$loss_raw_fro, 0.95, type = 7)),
      psd_fro_mean = psd_mean[["mean"]], psd_fro_mean_mcse = psd_mean[["mcse"]],
      psd_fro_rms = psd_rms[["rms"]], psd_fro_rms_mcse = psd_rms[["mcse"]],
      psd_fro_q05 = unname(quantile(x$loss_psd_fro, 0.05, type = 7)),
      psd_fro_median = median(x$loss_psd_fro),
      psd_fro_q95 = unname(quantile(x$loss_psd_fro, 0.95, type = 7)),
      stable_fro_mean = stable_mean[["mean"]], stable_fro_mean_mcse = stable_mean[["mcse"]],
      stable_fro_rms = stable_rms[["rms"]], stable_fro_rms_mcse = stable_rms[["mcse"]],
      loss_ratio_mean = ratio[["mean"]], loss_ratio_mean_mcse = ratio[["mcse"]],
      loss_ratio_max = max(x$loss_ratio_psd_raw),
      projection_improves_frequency = mean(x$projection_improves_loss),
      theorem_violations = sum(!x$projection_inequality_pass),
      minimum_theorem_margin = min(x$projection_inequality_margin),
      R_max_error_mean = R_max[["mean"]], R_max_error_mean_mcse = R_max[["mcse"]],
      R_fro_error_mean = R_fro[["mean"]], R_fro_error_mean_mcse = R_fro[["mcse"]],
      R_rmse_offdiag_mean = R_rmse[["mean"]], R_rmse_offdiag_mean_mcse = R_rmse[["mcse"]],
      partial_max_error_mean = P_max[["mean"]], partial_max_error_mean_mcse = P_max[["mcse"]],
      partial_fro_error_mean = P_fro[["mean"]], partial_fro_error_mean_mcse = P_fro[["mcse"]],
      partial_rmse_offdiag_mean = P_rmse[["mean"]],
      partial_rmse_offdiag_mean_mcse = P_rmse[["mcse"]],
      partial_max_error_gt_0_5 = mean(x$max_abs_partial_error > 0.5),
      partial_max_error_gt_0_75 = mean(x$max_abs_partial_error > 0.75),
      partial_rmse_mean_active = if (any(active)) mean(x$rmse_offdiag_partial[active]) else NA_real_,
      partial_rmse_mean_inactive = if (any(!active)) mean(x$rmse_offdiag_partial[!active]) else NA_real_,
      partial_max_mean_projection_active = if (any(x$projection_active)) mean(x$max_abs_partial_error[x$projection_active]) else NA_real_,
      partial_max_mean_projection_inactive = if (any(!x$projection_active)) mean(x$max_abs_partial_error[!x$projection_active]) else NA_real_,
      partial_max_mean_floor_active = if (any(x$floor_active)) mean(x$max_abs_partial_error[x$floor_active]) else NA_real_,
      partial_max_mean_floor_inactive = if (any(!x$floor_active)) mean(x$max_abs_partial_error[!x$floor_active]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

summarise_matrix_pairs <- function(results, type = c("correlation", "partial")) {
  type <- match.arg(type)
  validate_matrix_full_output(results)
  truth <- matrix_truth()
  pairs <- which(upper.tri(truth$R), arr.ind = TRUE)
  rows <- list()
  z <- 0L
  for (n in sort(unique(results$n))) {
    x <- results[results$n == n, , drop = FALSE]
    estimate_column <- if (type == "correlation") x$R_estimate else x$Partial_estimate
    truth_matrix <- if (type == "correlation") truth$R else truth$Partial
    for (i in seq_len(nrow(pairs))) {
      z <- z + 1L
      j <- pairs[i, 1L]
      k <- pairs[i, 2L]
      estimate <- vapply(estimate_column, function(A) A[j, k], numeric(1))
      error <- estimate - truth_matrix[j, k]
      rmse <- sqrt(mean(error^2))
      rows[[z]] <- data.frame(
        n = n, variable_j = paste0("V", j), variable_k = paste0("V", k),
        truth = truth_matrix[j, k], mean_estimate = mean(estimate), bias = mean(error),
        mcse_bias = sd(error) / sqrt(length(error)), rmse = rmse,
        mcse_rmse = if (rmse == 0) 0 else sd(error^2) / (2 * rmse * sqrt(length(error))),
        mean_absolute_error = mean(abs(error)), median_absolute_error = median(abs(error)),
        abs_error_gt_0_5 = mean(abs(error) > 0.5), abs_error_gt_0_75 = mean(abs(error) > 0.75),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
