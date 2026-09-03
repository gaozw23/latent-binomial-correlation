.resolve_coordinate <- function(fit, coordinate, argument) {
  if (is.character(coordinate)) {
    if (length(coordinate) != 1L || is.na(coordinate) ||
        !coordinate %in% fit$coordinate_names) {
      stop(sprintf("%s must be one exact coordinate name", argument), call. = FALSE)
    }
    return(match(coordinate, fit$coordinate_names))
  }
  if (!is.numeric(coordinate) || length(coordinate) != 1L || is.na(coordinate) ||
      coordinate != as.integer(coordinate) || coordinate < 1L || coordinate > fit$d) {
    stop(sprintf("%s must be one integer index in 1:%d or an exact coordinate name",
                 argument, fit$d), call. = FALSE)
  }
  as.integer(coordinate)
}

.validate_ci_inputs <- function(fit, j, k, level, scale) {
  if (!inherits(fit, "latent_binomial_cov")) {
    stop("fit must inherit from 'latent_binomial_cov'", call. = FALSE)
  }
  if (is.null(fit$Psi_hat)) stop("fit$Psi_hat is required; refit with keep_influence=TRUE", call. = FALSE)
  j <- .resolve_coordinate(fit, j, "j")
  k <- .resolve_coordinate(fit, k, "k")
  if (j == k) stop("j and k must identify distinct coordinates", call. = FALSE)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be strictly between 0 and 1", call. = FALSE)
  }
  scale <- match.arg(scale, c("fisher", "raw"))
  list(j = j, k = k, level = level, scale = scale)
}

.make_ci_result <- function(estimate, phi_raw, n, level, scale, diagnostics) {
  phi_raw <- phi_raw - mean(phi_raw)
  se_raw <- sqrt(mean(phi_raw^2) / n)
  if (abs(estimate) >= 1) {
    stop("Fisher inference requires an estimate strictly between -1 and 1", call. = FALSE)
  }
  z_estimate <- atanh(estimate)
  phi_fisher <- phi_raw / (1 - estimate^2)
  phi_fisher <- phi_fisher - mean(phi_fisher)
  se_fisher <- sqrt(mean(phi_fisher^2) / n)
  critical <- qnorm(1 - (1 - level) / 2)
  if (identical(scale, "fisher")) {
    limits <- tanh(z_estimate + c(-1, 1) * critical * se_fisher)
  } else {
    limits <- estimate + c(-1, 1) * critical * se_raw
  }
  list(
    estimate = estimate,
    z_estimate = z_estimate,
    se_raw = se_raw,
    se_fisher = se_fisher,
    lower = limits[1L],
    upper = limits[2L],
    level = level,
    scale = scale,
    phi_raw = phi_raw,
    phi_fisher = phi_fisher,
    regularity_diagnostics = diagnostics
  )
}

latent_correlation_ci <- function(fit, j, k, level = 0.95,
                                  scale = "fisher") {
  input <- .validate_ci_inputs(fit, j, k, level, scale)
  j <- input$j
  k <- input$k
  vj <- fit$Sigma_stable[j, j]
  vk <- fit$Sigma_stable[k, k]
  if (!is.finite(vj) || !is.finite(vk) || vj <= 0 || vk <= 0) {
    stop("Stable variances must be finite and strictly positive", call. = FALSE)
  }
  rho <- fit$R_latent[j, k]
  if (fit$projection_active) {
    warning("PSD projection is active; regular influence-function inference may be unreliable",
            call. = FALSE)
  }
  if (abs(rho) > 0.95) {
    warning("Absolute correlation exceeds 0.95; Fisher-scale regularity is near the boundary",
            call. = FALSE)
  }

  phi <- fit$Psi_hat[, j, k] / sqrt(vj * vk) -
    0.5 * rho * (fit$Psi_hat[, j, j] / vj + fit$Psi_hat[, k, k] / vk)
  diagnostics <- list(
    coordinate_j = fit$coordinate_names[j],
    coordinate_k = fit$coordinate_names[k],
    stable_variances_positive = TRUE,
    projection_active = fit$projection_active,
    floor_active = fit$floor_active,
    near_boundary = abs(rho) > 0.95,
    regular_inference_requires_positive_definite_interior_truth = TRUE
  )
  .make_ci_result(rho, phi, fit$n, input$level, input$scale, diagnostics)
}
