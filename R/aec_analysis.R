.aec_required_columns <- c(
  "match_id", "match_route", "State", "PollingPlaceID_2019",
  "PollingPlaceID_2022", "X_2019", "M_2019", "Y_2019",
  "X_2022", "M_2022", "Y_2022"
)

read_frozen_aec <- function(path = "results/derived/aec_matched_frozen.csv") {
  if (!file.exists(path)) stop("Frozen AEC match table is missing", call. = FALSE)
  out <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(.aec_required_columns, names(out))
  if (length(missing)) stop("Frozen AEC table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  numeric_columns <- c("X_2019", "M_2019", "Y_2019", "X_2022", "M_2022", "Y_2022")
  out[numeric_columns] <- lapply(out[numeric_columns], as.numeric)
  out
}

.aec_matrices <- function(data) {
  X <- cbind(`2019` = data$X_2019, `2022` = data$X_2022)
  M <- cbind(`2019` = data$M_2019, `2022` = data$M_2022)
  list(X = X, M = M, Y = X / M)
}

.aec_distribution <- function(x, prefix) {
  probabilities <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
  values <- c(stats::quantile(x, probabilities, type = 7, names = FALSE), mean(x))
  names(values) <- paste0(prefix, "_", c("min", "q05", "q25", "median", "q75", "q95", "max", "mean"))
  values
}

.aec_matrix_cells <- function(A, prefix) {
  values <- c(A[1, 1], A[1, 2], A[2, 1], A[2, 2])
  names(values) <- paste0(prefix, c("_11", "_12", "_21", "_22"))
  values
}

.aec_associations <- function(data) {
  values <- c(
    cor(data$M_2019, data$Y_2019, method = "pearson"),
    cor(data$M_2019, data$Y_2022, method = "pearson"),
    cor(data$M_2022, data$Y_2019, method = "pearson"),
    cor(data$M_2022, data$Y_2022, method = "pearson"),
    cor(data$M_2019, data$Y_2019, method = "spearman"),
    cor(data$M_2019, data$Y_2022, method = "spearman"),
    cor(data$M_2022, data$Y_2019, method = "spearman"),
    cor(data$M_2022, data$Y_2022, method = "spearman")
  )
  names(values) <- c(
    "pearson_M2019_Y2019", "pearson_M2019_Y2022",
    "pearson_M2022_Y2019", "pearson_M2022_Y2022",
    "spearman_M2019_Y2019", "spearman_M2019_Y2022",
    "spearman_M2022_Y2019", "spearman_M2022_Y2022"
  )
  values
}

analyse_aec <- function(data, analysis_set = "primary", primary_latent = NA_real_,
                        small_place_threshold = NA_real_) {
  matrices <- .aec_matrices(data)
  fit <- fit_latent_binomial_cov(matrices$X, matrices$M, keep_influence = TRUE)
  interval <- latent_correlation_ci(fit, 1L, 2L, level = 0.95, scale = "fisher")
  naive <- cor(matrices$Y[, 1], matrices$Y[, 2])
  latent <- fit$R_latent[1, 2]
  attenuation_factor <- sqrt(prod(fit$reliability))
  state_counts <- table(data$State)
  route_counts <- table(data$match_route)

  values <- c(
    analysis_set = analysis_set,
    n = nrow(data),
    exact_id_n = sum(data$match_route == "EXACT_ID"),
    secondary_n = sum(data$match_route == "UNIQUE_NAME_COORD"),
    small_place_threshold = small_place_threshold,
    .aec_distribution(data$M_2019, "M_2019"),
    .aec_distribution(data$M_2022, "M_2022"),
    .aec_distribution(data$Y_2019, "Y_2019"),
    .aec_distribution(data$Y_2022, "Y_2022"),
    naive_correlation = naive,
    .aec_matrix_cells(fit$Sigma_raw, "Sigma_raw"),
    .aec_matrix_cells(fit$Sigma_psd, "Sigma_psd"),
    .aec_matrix_cells(fit$Sigma_stable, "Sigma_stable"),
    min_eigen_raw = fit$min_eigen_raw,
    lambda_n = fit$lambda_n,
    projection_active = fit$projection_active,
    floor_active = fit$floor_active,
    latent_correlation = latent,
    reliability_2019 = unname(fit$reliability[1]),
    reliability_2022 = unname(fit$reliability[2]),
    attenuation_factor = attenuation_factor,
    fitted_attenuated_correlation = latent * attenuation_factor,
    analytic_se_fisher = interval$se_fisher,
    analytic_ci_lower = interval$lower,
    analytic_ci_upper = interval$upper,
    latent_minus_naive = latent - naive,
    relative_attenuation = if (abs(latent) > sqrt(.Machine$double.eps)) (latent - naive) / abs(latent) else NA_real_,
    difference_from_primary = if (is.finite(primary_latent)) latent - primary_latent else 0,
    state_counts = paste(names(state_counts), as.integer(state_counts), sep = "=", collapse = "|"),
    match_route_counts = paste(names(route_counts), as.integer(route_counts), sep = "=", collapse = "|"),
    .aec_associations(data)
  )
  result <- as.data.frame(as.list(values), stringsAsFactors = FALSE, check.names = FALSE)
  logical_names <- c("projection_active", "floor_active")
  result[logical_names] <- lapply(result[logical_names], as.logical)
  numeric_names <- setdiff(names(result), c(
    "analysis_set", "state_counts", "match_route_counts", logical_names
  ))
  result[numeric_names] <- lapply(result[numeric_names], as.numeric)
  list(summary = result, fit = fit, interval = interval)
}

run_aec_state_bootstrap <- function(data, streams) {
  if (length(streams) != 4999L) stop("AEC bootstrap requires exactly 4999 streams", call. = FALSE)
  matrices <- .aec_matrices(data)
  state_indices <- split(seq_len(nrow(data)), data$State)
  small_states <- names(state_indices)[lengths(state_indices) < 2L]
  one <- function(b) {
    with_rng_stream(streams[[b]], {
      index <- unlist(lapply(state_indices, function(z) {
        if (length(z) < 2L) z else sample(z, length(z), replace = TRUE)
      }), use.names = FALSE)
      result <- tryCatch({
        fit <- fit_latent_binomial_cov(
          matrices$X[index, , drop = FALSE], matrices$M[index, , drop = FALSE],
          keep_influence = FALSE
        )
        rho <- fit$R_latent[1, 2]
        valid <- is.finite(rho) && abs(rho) < 1
        list(rho = rho, fisher = if (valid) atanh(rho) else NA_real_,
             projection = fit$projection_active, floor = fit$floor_active,
             valid = valid, reason = if (valid) "" else "NONFINITE_OR_BOUNDARY")
      }, error = function(e) list(
        rho = NA_real_, fisher = NA_real_, projection = NA, floor = NA,
        valid = FALSE, reason = paste0("FIT_ERROR: ", conditionMessage(e))
      ))
      data.frame(
        replicate = b, rho = result$rho, fisher = result$fisher,
        projection_active = result$projection, floor_active = result$floor,
        valid = result$valid, failure_reason = result$reason,
        stringsAsFactors = FALSE
      )
    })
  }
  replicates <- do.call(rbind, lapply(seq_along(streams), one))
  valid <- replicates$valid & is.finite(replicates$fisher)
  z_limits <- quantile(replicates$fisher[valid], c(0.025, 0.975), type = 8, names = FALSE)
  rho_quantiles <- quantile(replicates$rho[valid], c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), type = 8, names = FALSE)
  summary <- data.frame(
    bootstrap_replicates = nrow(replicates), valid_replicates = sum(valid),
    invalid_replicates = sum(!valid), quantile_type = 8L,
    rho_min = rho_quantiles[1], rho_q025 = rho_quantiles[2], rho_q25 = rho_quantiles[3],
    rho_median = rho_quantiles[4], rho_mean = mean(replicates$rho[valid]),
    rho_q75 = rho_quantiles[5], rho_q975 = rho_quantiles[6], rho_max = rho_quantiles[7],
    bootstrap_ci_lower = tanh(z_limits[1]), bootstrap_ci_upper = tanh(z_limits[2]),
    projection_count = sum(replicates$projection_active %in% TRUE),
    floor_count = sum(replicates$floor_active %in% TRUE),
    deterministic_small_states = if (length(small_states)) paste(small_states, collapse = "|") else "NONE",
    stringsAsFactors = FALSE
  )
  list(replicates = replicates, summary = summary)
}

controlled_binomial_thinning <- function(data, streams,
                                         caps = c(5L, 10L, 20L, 50L, 100L),
                                         replicates_per_cap = 2000L) {
  expected <- length(caps) * replicates_per_cap
  if (length(streams) != expected) stop("Incorrect number of AEC thinning streams", call. = FALSE)
  pseudo_truth <- cor(data$Y_2019, data$Y_2022)
  mapping <- expand.grid(replicate = seq_len(replicates_per_cap), cap = caps,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  mapping <- mapping[order(mapping$cap, mapping$replicate), , drop = FALSE]
  mapping$stream_id <- seq_len(nrow(mapping))
  one <- function(i) {
    task <- mapping[i, ]
    with_rng_stream(streams[[task$stream_id]], {
      m19 <- pmin(data$M_2019, task$cap)
      m22 <- pmin(data$M_2022, task$cap)
      x19 <- rbinom(nrow(data), m19, data$Y_2019)
      x22 <- rbinom(nrow(data), m22, data$Y_2022)
      naive <- cor(x19 / m19, x22 / m22)
      result <- tryCatch({
        fit <- fit_latent_binomial_cov(cbind(`2019` = x19, `2022` = x22),
                                       cbind(`2019` = m19, `2022` = m22),
                                       keep_influence = FALSE)
        corrected <- fit$R_latent[1, 2]
        list(corrected = corrected, projection = fit$projection_active,
             floor = fit$floor_active, raw_invalid = any(diag(fit$Sigma_raw) <= 0),
             valid = is.finite(naive) && is.finite(corrected), reason = "")
      }, error = function(e) list(
        corrected = NA_real_, projection = NA, floor = NA, raw_invalid = NA,
        valid = FALSE, reason = paste0("FIT_ERROR: ", conditionMessage(e))
      ))
      if (!is.finite(naive)) {
        result$valid <- FALSE
        result$reason <- paste(result$reason, "NONFINITE_NAIVE")
      }
      data.frame(
        cap = task$cap, replicate = task$replicate, stream_id = task$stream_id,
        reference_correlation = pseudo_truth, naive = naive,
        corrected = result$corrected, projection_active = result$projection,
        floor_active = result$floor, invalid_raw_variance = result$raw_invalid,
        valid = result$valid, failure_reason = trimws(result$reason),
        stringsAsFactors = FALSE
      )
    })
  }
  results <- do.call(rbind, lapply(seq_len(nrow(mapping)), one))
  summaries <- lapply(caps, function(cap) {
    x <- results[results$cap == cap, , drop = FALSE]
    valid <- x$valid & is.finite(x$naive) & is.finite(x$corrected)
    naive_error <- x$naive[valid] - pseudo_truth
    corrected_error <- x$corrected[valid] - pseudo_truth
    data.frame(
      label = "SEMI-SYNTHETIC CONTROLLED THINNING", cap = cap,
      replicates = nrow(x), valid_replicates = sum(valid),
      failure_nonfinite_count = sum(!valid), reference_full_data_Y_correlation = pseudo_truth,
      mean_naive = mean(x$naive[valid]), naive_bias = mean(naive_error),
      naive_empirical_sd = sd(x$naive[valid]), naive_rmse = sqrt(mean(naive_error^2)),
      mean_corrected = mean(x$corrected[valid]), corrected_bias = mean(corrected_error),
      corrected_empirical_sd = sd(x$corrected[valid]), corrected_rmse = sqrt(mean(corrected_error^2)),
      projection_count = sum(x$projection_active %in% TRUE),
      projection_rate = mean(x$projection_active[valid]),
      floor_count = sum(x$floor_active %in% TRUE), floor_rate = mean(x$floor_active[valid]),
      invalid_raw_variance_count = sum(x$invalid_raw_variance %in% TRUE),
      stringsAsFactors = FALSE
    )
  })
  list(replicates = results, summary = do.call(rbind, summaries), mapping = mapping)
}
