testthat::test_that("T8 partial-correlation influence sign matches numerical differentiation", {
  Y <- rbind(c(0.10, 0.20, 0.30), c(0.20, 0.70, 0.40), c(0.50, 0.30, 0.80),
             c(0.80, 0.60, 0.20), c(0.70, 0.90, 0.70), c(0.40, 0.50, 0.10),
             c(0.30, 0.15, 0.65), c(0.90, 0.40, 0.55))
  Q <- matrix(c(
    0.002,0.003,0.002, 0.003,0.002,0.003, 0.002,0.002,0.003, 0.003,0.003,0.002,
    0.002,0.003,0.003, 0.003,0.002,0.002, 0.002,0.003,0.002, 0.003,0.002,0.003
  ), nrow = 8, byrow = TRUE)
  weights <- c(0.08, 0.12, 0.15, 0.16, 0.14, 0.13, 0.10, 0.12)
  support_index <- 5L
  base <- weighted_covariance_functional(Y, Q, weights)
  Sigma <- base$Sigma
  testthat::expect_gt(min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values), 0)
  Theta <- solve(Sigma)
  pi12 <- -Theta[1, 2] / sqrt(Theta[1, 1] * Theta[2, 2])
  psi <- tcrossprod(Y[support_index, ] - base$mu) - diag(Q[support_index, ], 3L) - Sigma
  K <- -Theta %*% psi %*% Theta
  analytic <- -K[1, 2] / sqrt(Theta[1, 1] * Theta[2, 2]) -
    0.5 * pi12 * (K[1, 1] / Theta[1, 1] + K[2, 2] / Theta[2, 2])

  epsilons <- c(1e-5, 5e-6, 1e-6)
  numerical <- vapply(epsilons, function(epsilon) {
    perturbed_weights <- weights
    perturbed_weights[support_index] <- perturbed_weights[support_index] + epsilon
    perturbed <- weighted_covariance_functional(Y, Q, perturbed_weights)$Sigma
    theta_perturbed <- solve(perturbed)
    pi_perturbed <- -theta_perturbed[1, 2] /
      sqrt(theta_perturbed[1, 1] * theta_perturbed[2, 2])
    (pi_perturbed - pi12) / epsilon
  }, numeric(1))
  relative_error <- abs(numerical - analytic) / max(abs(analytic), 1e-15)
  testthat::expect_lt(relative_error[3L], 2e-4)
  testthat::expect_lt(abs(numerical[3L] - analytic), abs(numerical[1L] - analytic))
  cat(sprintf("T8 analytic derivative: %.17g; numerical: %.17g; relative error: %.17g\n",
              analytic, numerical[3L], relative_error[3L]))
  record_phase1_diagnostic("T8", list(epsilons = epsilons, analytic = analytic,
                                       numerical = numerical, relative_error = relative_error,
                                       tolerance = 2e-4))
})
