testthat::test_that("T2 hand-computed matrix fixture matches every definition", {
  X <- matrix(c(0, 1, 2, 3, 1, 2, 1, 0), nrow = 4, ncol = 2,
              dimnames = list(NULL, c("A", "B")))
  M <- matrix(c(2, 2, 4, 4, 2, 4, 2, 2), nrow = 4, ncol = 2,
              dimnames = list(NULL, c("A", "B")))
  expected_Y <- matrix(c(
    0.0000000000000000, 0.5000000000000000, 0.5000000000000000, 0.7500000000000000,
    0.5000000000000000, 0.5000000000000000, 0.5000000000000000, 0.0000000000000000
  ), nrow = 4, ncol = 2, dimnames = list(NULL, c("A", "B")))
  expected_Q <- matrix(c(
    0.0000000000000000, 0.2500000000000000, 0.08333333333333333, 0.0625000000000000,
    0.2500000000000000, 0.08333333333333333, 0.2500000000000000, 0.0000000000000000
  ), nrow = 4, ncol = 2, dimnames = list(NULL, c("A", "B")))
  expected_S <- matrix(c(
    0.09895833333333333, -0.05208333333333333,
    -0.05208333333333333, 0.06250000000000000
  ), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  expected_raw <- matrix(c(
    0.00000000000000000, -0.05208333333333333,
    -0.05208333333333333, -0.08333333333333333
  ), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  expected_emp <- matrix(c(
    -0.02473958333333333, -0.03906250000000000,
    -0.03906250000000000, -0.09895833333333333
  ), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))

  fit <- fit_latent_binomial_cov(X, M)
  testthat::expect_equal(fit$Y, expected_Y, tolerance = 1e-13)
  testthat::expect_equal(fit$Q, expected_Q, tolerance = 1e-13)
  testthat::expect_equal(fit$S_Y, expected_S, tolerance = 1e-13)
  testthat::expect_equal(fit$Sigma_raw, expected_raw, tolerance = 1e-13)
  testthat::expect_equal(fit$Sigma_emp_n, expected_emp, tolerance = 1e-13)
  max_error <- max(abs(c(fit$Y - expected_Y, fit$Q - expected_Q,
                         fit$S_Y - expected_S, fit$Sigma_raw - expected_raw,
                         fit$Sigma_emp_n - expected_emp)))
  cat(sprintf("T2 maximum fixture error: %.17g\n", max_error))
  record_phase1_diagnostic("T2", list(max_absolute_error = max_error, tolerance = 1e-13))
})

testthat::test_that("T9 row permutations leave aggregates and intervals invariant", {
  result <- with_rng_stream(core_test_streams[[9L]], {
    n <- 500L
    X <- matrix(sample.int(101L, n * 2L, replace = TRUE) - 1L, n, 2L,
                dimnames = list(NULL, c("V1", "V2")))
    M <- matrix(100L, n, 2L, dimnames = dimnames(X))
    fit1 <- fit_latent_binomial_cov(X, M)
    order <- sample.int(n)
    fit2 <- fit_latent_binomial_cov(X[order, , drop = FALSE], M[order, , drop = FALSE])
    ci1 <- latent_correlation_ci(fit1, 1, 2)
    ci2 <- latent_correlation_ci(fit2, 1, 2)
    list(fit1 = fit1, fit2 = fit2, ci1 = ci1, ci2 = ci2, order = order)
  })
  aggregate_names <- c("ybar", "qbar", "S_Y", "Sigma_emp_n", "Sigma_raw",
                       "Sigma_psd", "Sigma_stable", "R_latent", "Theta", "Partial",
                       "reliability")
  errors <- vapply(aggregate_names, function(name) {
    max(abs(result$fit1[[name]] - result$fit2[[name]]))
  }, numeric(1))
  ci_error <- max(abs(unlist(result$ci1[c("estimate", "z_estimate", "se_raw", "se_fisher", "lower", "upper")]) -
                      unlist(result$ci2[c("estimate", "z_estimate", "se_raw", "se_fisher", "lower", "upper")])))
  testthat::expect_lte(max(errors), 1e-13)
  testthat::expect_lte(ci_error, 1e-13)
  cat(sprintf("T9 max aggregate error: %.17g; max CI error: %.17g\n", max(errors), ci_error))
  record_phase1_diagnostic("T9", list(max_aggregate_error = max(errors), max_ci_error = ci_error,
                                       tolerance = 1e-13))
})

testthat::test_that("T10 coordinate permutations are equivariant", {
  result <- with_rng_stream(core_test_streams[[10L]], {
    n <- 600L
    X <- matrix(sample.int(101L, n * 5L, replace = TRUE) - 1L, n, 5L,
                dimnames = list(NULL, paste0("V", seq_len(5L))))
    M <- matrix(100L, n, 5L, dimnames = dimnames(X))
    fit1 <- fit_latent_binomial_cov(X, M)
    perm <- c(3L, 1L, 5L, 2L, 4L)
    fit2 <- fit_latent_binomial_cov(X[, perm, drop = FALSE], M[, perm, drop = FALSE])
    list(fit1 = fit1, fit2 = fit2, perm = perm)
  })
  matrix_names <- c("S_Y", "Sigma_emp_n", "Sigma_raw", "Sigma_psd", "Sigma_stable",
                    "R_latent", "Theta", "Partial")
  matrix_errors <- vapply(matrix_names, function(name) {
    max(abs(unname(result$fit2[[name]]) - unname(result$fit1[[name]][result$perm, result$perm])))
  }, numeric(1))
  psi_expected <- result$fit1$Psi_hat[, result$perm, result$perm, drop = FALSE]
  psi_error <- max(abs(unname(result$fit2$Psi_hat) - unname(psi_expected)))
  vector_error <- max(abs(unname(result$fit2$reliability) - unname(result$fit1$reliability[result$perm])))
  testthat::expect_lte(max(matrix_errors), 1e-13)
  testthat::expect_lte(psi_error, 1e-13)
  testthat::expect_lte(vector_error, 1e-13)
  cat(sprintf("T10 max matrix error: %.17g; influence error: %.17g\n",
              max(matrix_errors), psi_error))
  record_phase1_diagnostic("T10", list(max_matrix_error = max(matrix_errors),
                                        max_influence_error = psi_error,
                                        max_vector_error = vector_error,
                                        tolerance = 1e-13))
})

testthat::test_that("T11 uses the exact n^-2 eigenvalue floor and not ridge addition", {
  n_values <- c(4L, 10L, 50L, 200L)
  floor_errors <- numeric(length(n_values))
  reconstruction_errors <- numeric(length(n_values))
  ridge_differences <- numeric(length(n_values))
  scaled <- numeric(length(n_values))
  for (ii in seq_along(n_values)) {
    n <- n_values[ii]
    i <- seq_len(n)
    M <- matrix(10, n, 3, dimnames = list(NULL, c("A", "B", "C")))
    X <- cbind(A = i %% 11, B = (3 * i + 1) %% 11, C = (7 * i + 2) %% 11)
    fit <- fit_latent_binomial_cov(X, M)
    exact_floor <- n^(-2)
    floor_errors[ii] <- abs(fit$lambda_n - exact_floor)
    e <- eigen(symmetrise_matrix(fit$Sigma_raw), symmetric = TRUE)
    expected <- e$vectors %*% diag(pmax(e$values, exact_floor), 3, 3) %*% t(e$vectors)
    expected <- symmetrise_matrix(expected)
    reconstruction_errors[ii] <- norm(unname(fit$Sigma_stable) - expected, "F")
    ridge_differences[ii] <- norm(unname(fit$Sigma_stable) -
                                    (unname(fit$Sigma_psd) + diag(exact_floor, 3)), "F")
    scaled[ii] <- sqrt(n) * exact_floor
    testthat::expect_identical(fit$lambda_n, exact_floor)
  }
  direct <- project_psd(diag(c(-0.3, 0.2)), eigen_floor = 0.01)
  testthat::expect_equal(direct$adjusted_eigenvalues, c(0.2, 0.01), tolerance = 0)
  testthat::expect_equal(sort(diag(direct$matrix)), c(0.01, 0.2), tolerance = 1e-15)
  testthat::expect_lte(max(reconstruction_errors), 1e-13)
  testthat::expect_true(any(ridge_differences > 1e-8))
  testthat::expect_true(all(diff(scaled) < 0))
  cat(sprintf("T11 max floor error: %.17g; max reconstruction error: %.17g\n",
              max(floor_errors), max(reconstruction_errors)))
  record_phase1_diagnostic("T11", list(n_values = n_values, floor_errors = floor_errors,
                                        scaled_sqrt_n_lambda = scaled,
                                        reconstruction_errors = reconstruction_errors,
                                        ridge_differences = ridge_differences))
})
