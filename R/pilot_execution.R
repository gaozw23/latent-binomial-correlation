.pilot_project_root <- function() {
  normalizePath(getOption("latent_binomial_project_root", "."), winslash = "/", mustWork = TRUE)
}

generate_phase2_stream_files <- function(root = .pilot_project_root()) {
  scenario <- read_pairwise_scenarios(file.path(root, "config", "pairwise_scenarios.csv"))
  pairwise_map <- do.call(rbind, lapply(scenario$scenario_id, function(id) {
    do.call(rbind, lapply(PAIRWISE_N_VALUES, function(n) {
      data.frame(scenario_id = id, n = n,
                 replicate = seq_len(PAIRWISE_PILOT_REPLICATES), stringsAsFactors = FALSE)
    }))
  }))
  pairwise_map$stream_id <- seq_len(nrow(pairwise_map))
  pairwise_map$study_component <- "pairwise_pilot"
  pairwise_map <- pairwise_map[c("stream_id", "scenario_id", "n", "replicate", "study_component")]
  pairwise_streams <- make_streams(nrow(pairwise_map), MASTER_SEEDS[["pairwise"]])

  matrix_map <- do.call(rbind, lapply(MATRIX_N_VALUES, function(n) {
    data.frame(n = n, replicate = seq_len(MATRIX_PILOT_REPLICATES), stringsAsFactors = FALSE)
  }))
  matrix_map$stream_id <- seq_len(nrow(matrix_map))
  matrix_map$study_component <- "matrix_pilot"
  matrix_map <- matrix_map[c("stream_id", "n", "replicate", "study_component")]
  matrix_streams <- make_streams(nrow(matrix_map), MASTER_SEEDS[["matrix"]])

  rmse_map <- expand.grid(
    scenario_id = scenario$scenario_id,
    n = PAIRWISE_N_VALUES,
    method = sort(c("proposed", "naive", "raw", "oracle_sample")),
    stringsAsFactors = FALSE
  )
  rmse_map <- rmse_map[order(rmse_map$scenario_id, rmse_map$n, rmse_map$method), , drop = FALSE]
  rownames(rmse_map) <- NULL
  rmse_map$stream_id <- seq_len(nrow(rmse_map))
  rmse_map$study_component <- "pilot_rmse_mcse"
  rmse_map <- rmse_map[c("stream_id", "scenario_id", "n", "method", "study_component")]
  rmse_streams <- make_streams(nrow(rmse_map), MASTER_SEEDS[["rmse_mcse"]])

  stream_dir <- file.path(root, "config", "rng_streams")
  dir.create(stream_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(
    pairwise_rds = file.path(stream_dir, "pairwise_pilot_streams.rds"),
    pairwise_csv = file.path(stream_dir, "pairwise_pilot_stream_map.csv"),
    matrix_rds = file.path(stream_dir, "matrix_pilot_streams.rds"),
    matrix_csv = file.path(stream_dir, "matrix_pilot_stream_map.csv"),
    rmse_rds = file.path(stream_dir, "rmse_pilot_streams.rds"),
    rmse_csv = file.path(stream_dir, "rmse_pilot_stream_map.csv")
  )
  atomic_save_rds(pairwise_streams, files[["pairwise_rds"]])
  write.csv(pairwise_map, files[["pairwise_csv"]], row.names = FALSE, quote = TRUE)
  atomic_save_rds(matrix_streams, files[["matrix_rds"]])
  write.csv(matrix_map, files[["matrix_csv"]], row.names = FALSE, quote = TRUE)
  atomic_save_rds(rmse_streams, files[["rmse_rds"]])
  write.csv(rmse_map, files[["rmse_csv"]], row.names = FALSE, quote = TRUE)

  stopifnot(length(readRDS(files[["pairwise_rds"]])) == 4800L,
            nrow(read.csv(files[["pairwise_csv"]])) == 4800L,
            length(readRDS(files[["matrix_rds"]])) == 1500L,
            nrow(read.csv(files[["matrix_csv"]])) == 1500L,
            length(readRDS(files[["rmse_rds"]])) == 96L,
            nrow(read.csv(files[["rmse_csv"]])) == 96L)
  data.frame(
    file = vapply(files, function(path) gsub("\\\\", "/", sub(paste0("^", root, "/?"), "", path)), character(1)),
    sha256 = vapply(files, sha256_file, character(1)),
    bytes = vapply(files, function(path) file.info(path)$size, numeric(1)),
    stringsAsFactors = FALSE
  )
}

atomic_save_rds <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  saveRDS(object, temporary, compress = compress, version = 3)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) stop("Could not atomically move RDS checkpoint into place", call. = FALSE)
  invisible(path)
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Package 'digest' is required", call. = FALSE)
  digest::digest(path, algo = "sha256", file = TRUE)
}

canonical_pairwise_tasks <- function(streams, mapping) {
  required <- c("stream_id", "scenario_id", "n", "replicate", "study_component")
  if (!all(required %in% names(mapping))) stop("Invalid pairwise RNG stream mapping schema", call. = FALSE)
  if (length(streams) != nrow(mapping)) stop("Pairwise stream count does not match mapping", call. = FALSE)
  order <- order(mapping$scenario_id, mapping$n, mapping$replicate)
  mapping <- mapping[order, , drop = FALSE]
  rownames(mapping) <- NULL
  lapply(seq_len(nrow(mapping)), function(i) {
    list(stream_id = mapping$stream_id[i], scenario_id = mapping$scenario_id[i],
         n = as.integer(mapping$n[i]), replicate = as.integer(mapping$replicate[i]),
         seed = streams[[mapping$stream_id[i]]])
  })
}

canonical_matrix_tasks <- function(streams, mapping) {
  required <- c("stream_id", "n", "replicate", "study_component")
  if (!all(required %in% names(mapping))) stop("Invalid matrix RNG stream mapping schema", call. = FALSE)
  if (length(streams) != nrow(mapping)) stop("Matrix stream count does not match mapping", call. = FALSE)
  order <- order(mapping$n, mapping$replicate)
  mapping <- mapping[order, , drop = FALSE]
  rownames(mapping) <- NULL
  lapply(seq_len(nrow(mapping)), function(i) {
    list(stream_id = mapping$stream_id[i], n = as.integer(mapping$n[i]),
         replicate = as.integer(mapping$replicate[i]), seed = streams[[mapping$stream_id[i]]])
  })
}

open_pilot_cluster <- function(workers) {
  workers <- as.integer(workers)
  if (length(workers) != 1L || is.na(workers) || workers < 1L) stop("workers must be a positive integer", call. = FALSE)
  if (workers == 1L) return(NULL)
  cluster <- parallel::makePSOCKcluster(workers)
  root <- .pilot_project_root()
  parallel::clusterExport(cluster, "root", envir = environment())
  setup <- parallel::clusterEvalQ(cluster, {
    options(latent_binomial_project_root = root)
    files <- sort(list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE))
    invisible(lapply(files, sys.source, envir = .GlobalEnv))
    TRUE
  })
  if (!all(unlist(setup))) {
    parallel::stopCluster(cluster)
    stop("Pilot worker initialization failed", call. = FALSE)
  }
  cluster
}

close_pilot_cluster <- function(cluster) {
  if (!is.null(cluster)) parallel::stopCluster(cluster)
  invisible(TRUE)
}

execute_pairwise_tasks <- function(tasks, cluster = NULL, load_balanced = TRUE) {
  if (!length(tasks)) return(data.frame())
  rows <- if (is.null(cluster)) {
    lapply(tasks, pairwise_pilot_replicate)
  } else if (isTRUE(load_balanced)) {
    parallel::parLapplyLB(cluster, tasks, function(task) pairwise_pilot_replicate(task))
  } else {
    parallel::parLapply(cluster, tasks, function(task) pairwise_pilot_replicate(task))
  }
  canonical_pairwise_results(do.call(rbind, rows))
}

execute_matrix_tasks <- function(tasks, cluster = NULL, load_balanced = TRUE) {
  if (!length(tasks)) return(data.frame())
  rows <- if (is.null(cluster)) {
    lapply(tasks, matrix_pilot_replicate)
  } else if (isTRUE(load_balanced)) {
    parallel::parLapplyLB(cluster, tasks, function(task) matrix_pilot_replicate(task))
  } else {
    parallel::parLapply(cluster, tasks, function(task) matrix_pilot_replicate(task))
  }
  canonical_matrix_results(do.call(rbind, rows))
}

.pairwise_result_template <- function(task) {
  data.frame(
    scenario_id = as.character(task$scenario_id), n = as.integer(task$n),
    replicate = as.integer(task$replicate), stream_id = as.integer(task$stream_id),
    rho_true = NA_real_, rho_proposed = NA_real_, rho_raw = NA_real_,
    rho_raw_valid = FALSE, rho_naive = NA_real_, rho_oracle_sample = NA_real_,
    if_lower = NA_real_, if_upper = NA_real_, if_cover = NA,
    if_length = NA_real_, if_se_fisher = NA_real_,
    Sigma11_raw = NA_real_, Sigma22_raw = NA_real_, Sigma12_raw = NA_real_,
    Sigma11_psd = NA_real_, Sigma22_psd = NA_real_, Sigma12_psd = NA_real_,
    min_eigen_raw = NA_real_, projection_active = NA, floor_active = NA,
    negative_raw_variance = NA, mean_M1 = NA_real_, mean_M2 = NA_real_,
    median_M1 = NA_real_, median_M2 = NA_real_, cor_M1_M2 = NA_real_,
    cor_M1_P1 = NA_real_, cor_M2_P2 = NA_real_,
    success = FALSE, warning_count = 0L, warning_messages = "",
    nonfinite_count = NA_integer_, failure_message = "",
    probability_bounds_valid = FALSE, denominator_valid = FALSE,
    count_valid = FALSE, denominator_mechanism_valid = FALSE,
    unweighted_valid = FALSE, rho_truth_analytic = FALSE,
    stringsAsFactors = FALSE
  )
}

pairwise_pilot_replicate <- function(task) {
  out <- .pairwise_result_template(task)
  warnings <- character()
  result <- tryCatch(
    withCallingHandlers({
      assign(".Random.seed", task$seed, envir = .GlobalEnv)
      dat <- generate_pairwise_data(task$n, task$scenario_id)
      fit <- fit_latent_binomial_cov(dat$X, dat$M, keep_influence = TRUE)
      ci <- latent_correlation_ci(fit, 1L, 2L, level = 0.95, scale = "fisher")

      raw_valid <- all(diag(fit$Sigma_raw) > 0)
      rho_raw <- if (raw_valid) {
        fit$Sigma_raw[1L, 2L] / sqrt(fit$Sigma_raw[1L, 1L] * fit$Sigma_raw[2L, 2L])
      } else NA_real_
      rho_naive <- cor(fit$Y[, 1L], fit$Y[, 2L])
      rho_oracle <- cor(dat$P[, 1L], dat$P[, 2L])

      probability_ok <- all(is.finite(dat$P)) && all(dat$P >= 0 & dat$P <= 1)
      denominator_ok <- all(dat$M >= 2) && all(abs(dat$M - round(dat$M)) <= INTEGER_TOL)
      count_ok <- all(abs(dat$X - round(dat$X)) <= INTEGER_TOL) &&
        all(dat$X >= 0 & dat$X <= dat$M)
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
      analytic_truth <- identical(unname(dat$truth$rho), unname(pairwise_truth(task$scenario_id)$rho))

      required_finite <- c(fit$R_latent[1L, 2L], rho_naive, rho_oracle,
                           ci$lower, ci$upper, ci$se_fisher,
                           fit$Sigma_raw, fit$Sigma_psd, fit$Sigma_stable,
                           fit$Theta, fit$Partial)
      nonfinite <- sum(!is.finite(required_finite))

      out$rho_true <- dat$truth$rho
      out$rho_proposed <- fit$R_latent[1L, 2L]
      out$rho_raw <- rho_raw
      out$rho_raw_valid <- raw_valid
      out$rho_naive <- rho_naive
      out$rho_oracle_sample <- rho_oracle
      out$if_lower <- ci$lower
      out$if_upper <- ci$upper
      out$if_cover <- ci$lower <= dat$truth$rho && ci$upper >= dat$truth$rho
      out$if_length <- ci$upper - ci$lower
      out$if_se_fisher <- ci$se_fisher
      out$Sigma11_raw <- fit$Sigma_raw[1L, 1L]
      out$Sigma22_raw <- fit$Sigma_raw[2L, 2L]
      out$Sigma12_raw <- fit$Sigma_raw[1L, 2L]
      out$Sigma11_psd <- fit$Sigma_psd[1L, 1L]
      out$Sigma22_psd <- fit$Sigma_psd[2L, 2L]
      out$Sigma12_psd <- fit$Sigma_psd[1L, 2L]
      out$min_eigen_raw <- fit$min_eigen_raw
      out$projection_active <- fit$projection_active
      out$floor_active <- fit$floor_active
      out$negative_raw_variance <- any(diag(fit$Sigma_raw) < 0)
      out$mean_M1 <- mean(dat$M[, 1L])
      out$mean_M2 <- mean(dat$M[, 2L])
      out$median_M1 <- median(dat$M[, 1L])
      out$median_M2 <- median(dat$M[, 2L])
      out$cor_M1_M2 <- cor(dat$M[, 1L], dat$M[, 2L])
      out$cor_M1_P1 <- cor(dat$M[, 1L], dat$P[, 1L])
      out$cor_M2_P2 <- cor(dat$M[, 2L], dat$P[, 2L])
      out$nonfinite_count <- as.integer(nonfinite)
      out$probability_bounds_valid <- probability_ok
      out$denominator_valid <- denominator_ok
      out$count_valid <- count_ok
      out$denominator_mechanism_valid <- mechanism_ok
      out$unweighted_valid <- isTRUE(fit$validation_summary$unweighted_rows)
      out$rho_truth_analytic <- analytic_truth
      out$success <- probability_ok && denominator_ok && count_ok && mechanism_ok &&
        analytic_truth && out$unweighted_valid && nonfinite == 0L
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

canonical_pairwise_results <- function(x) {
  if (!nrow(x)) return(x)
  x <- x[order(x$scenario_id, x$n, x$replicate), , drop = FALSE]
  rownames(x) <- NULL
  x
}

.matrix_result_template <- function(task) {
  data.frame(
    n = as.integer(task$n), replicate = as.integer(task$replicate),
    stream_id = as.integer(task$stream_id),
    loss_raw_fro = NA_real_, loss_psd_fro = NA_real_, loss_ratio_fro = NA_real_,
    loss_raw_operator = NA_real_, loss_psd_operator = NA_real_,
    min_eigen_raw = NA_real_, projection_active = NA, floor_active = NA,
    max_abs_R_error = NA_real_, rmse_offdiag_R = NA_real_,
    max_abs_partial_error = NA_real_, rmse_offdiag_partial = NA_real_,
    projection_inequality_pass = NA, projection_inequality_margin = NA_real_,
    success = FALSE, warning_count = 0L, warning_messages = "",
    nonfinite_count = NA_integer_, failure_message = "",
    probability_bounds_valid = FALSE, denominator_valid = FALSE,
    count_valid = FALSE, W_valid = FALSE, truth_valid = FALSE,
    stringsAsFactors = FALSE
  )
}

matrix_pilot_replicate <- function(task) {
  out <- .matrix_result_template(task)
  warnings <- character()
  result <- tryCatch(
    withCallingHandlers({
      assign(".Random.seed", task$seed, envir = .GlobalEnv)
      dat <- generate_matrix_data(task$n)
      fit <- fit_latent_binomial_cov(dat$X, dat$M, keep_influence = FALSE)
      truth <- dat$truth
      W_ok <- all(truth$W >= 0) && all(rowSums(truth$W) == 1)
      sigma_expected <- truth$sigma2_B * truth$W %*% t(truth$W)
      truth_ok <- isTRUE(all.equal(truth$Sigma, sigma_expected, tolerance = 0)) &&
        isTRUE(all.equal(truth$R, cov2cor(truth$Sigma), tolerance = 1e-15)) &&
        isTRUE(all.equal(truth$Theta, solve(truth$Sigma), tolerance = 1e-12))
      probability_ok <- all(is.finite(dat$P)) && all(dat$P >= 0 & dat$P <= 1)
      denominator_ok <- all(dat$M >= 2) && all(abs(dat$M - round(dat$M)) <= INTEGER_TOL)
      count_ok <- all(abs(dat$X - round(dat$X)) <= INTEGER_TOL) && all(dat$X >= 0 & dat$X <= dat$M)

      raw_error <- symmetrise_matrix(fit$Sigma_raw - truth$Sigma)
      psd_error <- symmetrise_matrix(fit$Sigma_psd - truth$Sigma)
      raw_fro <- norm(raw_error, "F")
      psd_fro <- norm(psd_error, "F")
      allowed <- raw_fro + PROJECTION_LOSS_TOL * (1 + raw_fro)
      margin <- allowed - psd_fro
      theorem_pass <- margin >= 0
      raw_operator <- max(abs(eigen(raw_error, symmetric = TRUE, only.values = TRUE)$values))
      psd_operator <- max(abs(eigen(psd_error, symmetric = TRUE, only.values = TRUE)$values))
      R_error <- fit$R_latent - truth$R
      partial_error <- fit$Partial - truth$Partial
      numeric_required <- c(raw_fro, psd_fro, raw_operator, psd_operator,
                            fit$min_eigen_raw, R_error, partial_error, margin,
                            fit$Sigma_raw, fit$Sigma_psd, fit$Sigma_stable)
      nonfinite <- sum(!is.finite(numeric_required))

      out$loss_raw_fro <- raw_fro
      out$loss_psd_fro <- psd_fro
      out$loss_ratio_fro <- if (raw_fro == 0) 0 else psd_fro / raw_fro
      out$loss_raw_operator <- raw_operator
      out$loss_psd_operator <- psd_operator
      out$min_eigen_raw <- fit$min_eigen_raw
      out$projection_active <- fit$projection_active
      out$floor_active <- fit$floor_active
      out$max_abs_R_error <- max(abs(R_error))
      out$rmse_offdiag_R <- sqrt(mean(R_error[upper.tri(R_error)]^2))
      out$max_abs_partial_error <- max(abs(partial_error))
      out$rmse_offdiag_partial <- sqrt(mean(partial_error[upper.tri(partial_error)]^2))
      out$projection_inequality_pass <- theorem_pass
      out$projection_inequality_margin <- margin
      out$nonfinite_count <- as.integer(nonfinite)
      out$probability_bounds_valid <- probability_ok
      out$denominator_valid <- denominator_ok
      out$count_valid <- count_ok
      out$W_valid <- W_ok
      out$truth_valid <- truth_ok
      out$success <- probability_ok && denominator_ok && count_ok && W_ok && truth_ok &&
        theorem_pass && nonfinite == 0L
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

canonical_matrix_results <- function(x) {
  if (!nrow(x)) return(x)
  x <- x[order(x$n, x$replicate), , drop = FALSE]
  rownames(x) <- NULL
  x
}

.task_key <- function(task, kind) {
  if (kind == "pairwise") paste(task$scenario_id, task$n, task$replicate, sep = "|")
  else paste(task$n, task$replicate, sep = "|")
}

.result_keys <- function(results, kind) {
  if (!nrow(results)) character()
  else if (kind == "pairwise") paste(results$scenario_id, results$n, results$replicate, sep = "|")
  else paste(results$n, results$replicate, sep = "|")
}

run_checkpointed_tasks <- function(tasks, kind = c("pairwise", "matrix"), checkpoint_path,
                                   final_path, cluster = NULL, checkpoint_every = 25L,
                                   max_new_tasks = Inf, load_balanced = TRUE) {
  kind <- match.arg(kind)
  existing <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else data.frame()
  completed <- .result_keys(existing, kind)
  pending <- tasks[!vapply(tasks, .task_key, character(1), kind = kind) %in% completed]
  if (is.finite(max_new_tasks)) pending <- head(pending, as.integer(max_new_tasks))
  if (length(pending)) {
    starts <- seq.int(1L, length(pending), by = as.integer(checkpoint_every))
    for (start in starts) {
      chunk <- pending[start:min(start + checkpoint_every - 1L, length(pending))]
      new <- if (kind == "pairwise") execute_pairwise_tasks(chunk, cluster, load_balanced)
      else execute_matrix_tasks(chunk, cluster, load_balanced)
      existing <- if (nrow(existing)) rbind(existing, new) else new
      existing <- if (kind == "pairwise") canonical_pairwise_results(existing)
      else canonical_matrix_results(existing)
      atomic_save_rds(existing, checkpoint_path)
    }
  }
  all_keys <- vapply(tasks, .task_key, character(1), kind = kind)
  complete <- setequal(.result_keys(existing, kind), all_keys) && nrow(existing) == length(tasks)
  if (complete) atomic_save_rds(existing, final_path)
  list(results = existing, complete = complete, completed = nrow(existing), planned = length(tasks),
       checkpoint_path = checkpoint_path, final_path = final_path)
}

pairwise_required_columns <- function() c(
  "scenario_id", "n", "replicate", "stream_id", "rho_true", "rho_proposed",
  "rho_raw", "rho_raw_valid", "rho_naive", "rho_oracle_sample", "if_lower",
  "if_upper", "if_cover", "if_length", "if_se_fisher", "Sigma11_raw",
  "Sigma22_raw", "Sigma12_raw", "Sigma11_psd", "Sigma22_psd", "Sigma12_psd",
  "min_eigen_raw", "projection_active", "floor_active", "negative_raw_variance",
  "mean_M1", "mean_M2", "median_M1", "median_M2", "cor_M1_M2", "cor_M1_P1",
  "cor_M2_P2"
)

matrix_required_columns <- function() c(
  "n", "replicate", "stream_id", "loss_raw_fro", "loss_psd_fro",
  "loss_ratio_fro", "loss_raw_operator", "loss_psd_operator", "min_eigen_raw",
  "projection_active", "floor_active", "max_abs_R_error", "rmse_offdiag_R",
  "max_abs_partial_error", "rmse_offdiag_partial", "projection_inequality_pass",
  "projection_inequality_margin"
)

validate_pairwise_pilot_output <- function(x, expected_rows = NULL) {
  missing <- setdiff(pairwise_required_columns(), names(x))
  if (length(missing)) stop("Pairwise output missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is.null(expected_rows) && nrow(x) != expected_rows) stop("Unexpected pairwise row count", call. = FALSE)
  keys <- paste(x$scenario_id, x$n, x$replicate, sep = "|")
  if (anyDuplicated(keys)) stop("Duplicate pairwise task keys", call. = FALSE)
  if (any(!x$success)) stop("Pairwise pilot contains failed tasks", call. = FALSE)
  if (sum(x$nonfinite_count) != 0L) stop("Pairwise pilot contains non-finite required results", call. = FALSE)
  generator_fields <- c("probability_bounds_valid", "denominator_valid", "count_valid",
                        "denominator_mechanism_valid", "unweighted_valid", "rho_truth_analytic")
  if (!all(vapply(x[generator_fields], all, logical(1)))) stop("Pairwise generator contract failed", call. = FALSE)
  invisible(TRUE)
}

validate_matrix_pilot_output <- function(x, expected_rows = NULL) {
  missing <- setdiff(matrix_required_columns(), names(x))
  if (length(missing)) stop("Matrix output missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is.null(expected_rows) && nrow(x) != expected_rows) stop("Unexpected matrix row count", call. = FALSE)
  keys <- paste(x$n, x$replicate, sep = "|")
  if (anyDuplicated(keys)) stop("Duplicate matrix task keys", call. = FALSE)
  if (any(!x$success)) stop("Matrix pilot contains failed tasks", call. = FALSE)
  if (sum(x$nonfinite_count) != 0L) stop("Matrix pilot contains non-finite required results", call. = FALSE)
  if (any(!x$projection_inequality_pass)) stop("Deterministic projection theorem violation", call. = FALSE)
  invisible(TRUE)
}

rmse_mcse_bootstrap <- function(errors, stream, resamples = 2000L) {
  errors <- errors[is.finite(errors)]
  if (!length(errors)) return(NA_real_)
  assign(".Random.seed", stream, envir = .GlobalEnv)
  values <- replicate(resamples, sqrt(mean(sample(errors^2, replace = TRUE))))
  sd(values)
}

summarise_pairwise_pilot <- function(results, rmse_streams, rmse_mapping) {
  validate_pairwise_pilot_output(results)
  methods <- c(proposed = "rho_proposed", naive = "rho_naive",
               raw = "rho_raw", oracle_sample = "rho_oracle_sample")
  groups <- unique(results[c("scenario_id", "n")])
  groups <- groups[order(groups$scenario_id, groups$n), , drop = FALSE]
  rows <- list()
  z <- 0L
  for (g in seq_len(nrow(groups))) {
    subset <- results[results$scenario_id == groups$scenario_id[g] & results$n == groups$n[g], , drop = FALSE]
    truth <- unique(subset$rho_true)
    if (length(truth) != 1L) stop("Non-unique analytic truth within pairwise cell", call. = FALSE)
    for (method in names(methods)) {
      z <- z + 1L
      estimate <- subset[[methods[[method]]]]
      valid <- is.finite(estimate)
      error <- estimate[valid] - truth
      stream_row <- rmse_mapping$scenario_id == groups$scenario_id[g] &
        rmse_mapping$n == groups$n[g] & rmse_mapping$method == method
      if (sum(stream_row) != 1L) stop("Missing RMSE MCSE stream mapping", call. = FALSE)
      stream <- rmse_streams[[rmse_mapping$stream_id[stream_row]]]
      rows[[z]] <- data.frame(
        scenario_id = groups$scenario_id[g], n = groups$n[g], method = method,
        B_total = nrow(subset), B_valid = sum(valid), truth = truth,
        mean_estimate = mean(estimate[valid]), bias = mean(error),
        mcse_bias = sd(error) / sqrt(length(error)), empirical_sd = sd(estimate[valid]),
        rmse = sqrt(mean(error^2)), mcse_rmse = rmse_mcse_bootstrap(error, stream),
        median_absolute_error = median(abs(error)),
        coverage = if (method == "proposed") mean(subset$if_cover) else NA_real_,
        mcse_coverage = if (method == "proposed") {
          p <- mean(subset$if_cover); sqrt(p * (1 - p) / nrow(subset))
        } else NA_real_,
        mean_interval_length = if (method == "proposed") mean(subset$if_length) else NA_real_,
        lower_noncoverage = if (method == "proposed") mean(subset$if_lower > truth) else NA_real_,
        upper_noncoverage = if (method == "proposed") mean(subset$if_upper < truth) else NA_real_,
        mean_estimated_se_fisher = if (method == "proposed") mean(subset$if_se_fisher) else NA_real_,
        projection_rate = mean(subset$projection_active),
        negative_variance_rate = mean(subset$negative_raw_variance),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarise_matrix_pilot <- function(results) {
  validate_matrix_pilot_output(results)
  rows <- lapply(sort(unique(results$n)), function(n) {
    x <- results[results$n == n, , drop = FALSE]
    data.frame(
      n = n, replicates = nrow(x), projection_rate = mean(x$projection_active),
      floor_rate = mean(x$floor_active), min_eigen_min = min(x$min_eigen_raw),
      min_eigen_q05 = unname(quantile(x$min_eigen_raw, 0.05, type = 7)),
      min_eigen_median = median(x$min_eigen_raw),
      min_eigen_q95 = unname(quantile(x$min_eigen_raw, 0.95, type = 7)),
      raw_fro_mean = mean(x$loss_raw_fro), psd_fro_mean = mean(x$loss_psd_fro),
      loss_ratio_mean = mean(x$loss_ratio_fro), loss_ratio_max = max(x$loss_ratio_fro),
      R_max_error_mean = mean(x$max_abs_R_error), R_rmse_mean = mean(x$rmse_offdiag_R),
      partial_max_error_mean = mean(x$max_abs_partial_error),
      partial_rmse_mean = mean(x$rmse_offdiag_partial),
      theorem_violations = sum(!x$projection_inequality_pass),
      minimum_theorem_margin = min(x$projection_inequality_margin),
      warning_count = sum(x$warning_count), nonfinite_count = sum(x$nonfinite_count),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_t12_reproducibility <- function(pairwise_tasks, output_dir = "results/audits") {
  tasks <- head(pairwise_tasks, 50L)
  one <- execute_pairwise_tasks(tasks, cluster = NULL, load_balanced = FALSE)
  cluster <- open_pilot_cluster(2L)
  on.exit(close_pilot_cluster(cluster), add = TRUE)
  two <- execute_pairwise_tasks(tasks, cluster = cluster, load_balanced = TRUE)
  validate_pairwise_pilot_output(one, 50L)
  validate_pairwise_pilot_output(two, 50L)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  one_path <- file.path(output_dir, "t12_workers1.rds")
  two_path <- file.path(output_dir, "t12_workers2.rds")
  atomic_save_rds(one, one_path)
  atomic_save_rds(two, two_path)
  list(pass = identical(one, two), workers_1_hash = sha256_file(one_path),
       workers_2_hash = sha256_file(two_path), tasks = 50L,
       load_balanced_parallel = TRUE)
}

run_checkpoint_resume_test <- function(cell_tasks, output_dir = "results/audits") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  uninterrupted_path <- file.path(output_dir, "checkpoint_uninterrupted.rds")
  checkpoint_path <- file.path(output_dir, "checkpoint_partial.rds")
  resumed_path <- file.path(output_dir, "checkpoint_resumed.rds")
  unlink(c(uninterrupted_path, checkpoint_path, resumed_path), force = TRUE)

  uninterrupted <- execute_pairwise_tasks(cell_tasks, cluster = NULL, load_balanced = FALSE)
  atomic_save_rds(uninterrupted, uninterrupted_path)

  cluster2 <- open_pilot_cluster(2L)
  interrupted <- run_checkpointed_tasks(
    cell_tasks, "pairwise", checkpoint_path, resumed_path, cluster2,
    checkpoint_every = 25L, max_new_tasks = 73L, load_balanced = TRUE
  )
  close_pilot_cluster(cluster2)
  if (interrupted$complete || interrupted$completed != 73L) {
    stop("Deliberate checkpoint interruption did not stop at 73 tasks", call. = FALSE)
  }

  cluster3 <- open_pilot_cluster(3L)
  resumed <- run_checkpointed_tasks(
    cell_tasks, "pairwise", checkpoint_path, resumed_path, cluster3,
    checkpoint_every = 25L, max_new_tasks = Inf, load_balanced = TRUE
  )
  close_pilot_cluster(cluster3)
  if (!resumed$complete) stop("Checkpoint resume did not complete the cell", call. = FALSE)
  validate_pairwise_pilot_output(resumed$results, length(cell_tasks))
  list(
    pass = identical(uninterrupted, resumed$results),
    tasks = length(cell_tasks), interrupted_after = 73L,
    uninterrupted_hash = sha256_file(uninterrupted_path),
    resumed_hash = sha256_file(resumed_path),
    checkpoint_hash_after_resume = sha256_file(checkpoint_path),
    resume_workers = 3L, interruption_workers = 2L
  )
}
