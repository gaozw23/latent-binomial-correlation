testthat::test_that("T7 correlation influence formula matches numerical differentiation", {
  Y <- rbind(c(0.10, 0.20), c(0.20, 0.70), c(0.50, 0.30),
             c(0.80, 0.60), c(0.70, 0.90), c(0.40, 0.50))
  Q <- rbind(c(0.004, 0.006), c(0.005, 0.004), c(0.003, 0.005),
             c(0.006, 0.003), c(0.004, 0.004), c(0.005, 0.006))
  weights <- c(0.10, 0.15, 0.20, 0.20, 0.15, 0.20)
  support_index <- 4L
  base <- weighted_covariance_functional(Y, Q, weights)
  Sigma <- base$Sigma
  rho <- Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2])
  psi <- tcrossprod(Y[support_index, ] - base$mu) - diag(Q[support_index, ], 2L) - Sigma
  analytic <- psi[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2]) -
    0.5 * rho * (psi[1, 1] / Sigma[1, 1] + psi[2, 2] / Sigma[2, 2])

  epsilons <- c(1e-5, 5e-6, 1e-6)
  numerical <- vapply(epsilons, function(epsilon) {
    perturbed_weights <- weights
    perturbed_weights[support_index] <- perturbed_weights[support_index] + epsilon
    perturbed <- weighted_covariance_functional(Y, Q, perturbed_weights)$Sigma
    rho_perturbed <- perturbed[1, 2] / sqrt(perturbed[1, 1] * perturbed[2, 2])
    (rho_perturbed - rho) / epsilon
  }, numeric(1))
  relative_error <- abs(numerical - analytic) / max(abs(analytic), 1e-15)
  testthat::expect_lt(relative_error[3L], 2e-4)
  testthat::expect_lt(abs(numerical[3L] - analytic), abs(numerical[1L] - analytic))
  cat(sprintf("T7 analytic derivative: %.17g; numerical: %.17g; relative error: %.17g\n",
              analytic, numerical[3L], relative_error[3L]))
  record_phase1_diagnostic("T7", list(epsilons = epsilons, analytic = analytic,
                                       numerical = numerical, relative_error = relative_error,
                                       tolerance = 2e-4))
})
