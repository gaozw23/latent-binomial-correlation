latent_partial_ci <- function(fit, j, k, level = 0.95,
                              scale = "fisher") {
  input <- .validate_ci_inputs(fit, j, k, level, scale)
  j <- input$j
  k <- input$k
  Theta <- fit$Theta
  if (any(!is.finite(Theta)) || Theta[j, j] <= 0 || Theta[k, k] <= 0) {
    stop("Precision diagonal entries must be finite and strictly positive", call. = FALSE)
  }
  pi_jk <- fit$Partial[j, k]
  if (fit$projection_active) {
    warning("PSD projection is active; regular partial-correlation inference may be unreliable",
            call. = FALSE)
  }
  if (abs(pi_jk) > 0.95) {
    warning("Absolute partial correlation exceeds 0.95; Fisher-scale regularity is near the boundary",
            call. = FALSE)
  }

  phi <- numeric(fit$n)
  denom <- sqrt(Theta[j, j] * Theta[k, k])
  for (i in seq_len(fit$n)) {
    K_i <- -Theta %*% fit$Psi_hat[i, , ] %*% Theta
    phi[i] <- -K_i[j, k] / denom -
      0.5 * pi_jk * (K_i[j, j] / Theta[j, j] + K_i[k, k] / Theta[k, k])
  }
  diagnostics <- list(
    coordinate_j = fit$coordinate_names[j],
    coordinate_k = fit$coordinate_names[k],
    precision_diagonals_positive = TRUE,
    projection_active = fit$projection_active,
    floor_active = fit$floor_active,
    near_boundary = abs(pi_jk) > 0.95,
    regular_inference_requires_positive_definite_interior_truth = TRUE
  )
  .make_ci_result(pi_jk, phi, fit$n, input$level, input$scale, diagnostics)
}
