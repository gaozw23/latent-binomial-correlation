testthat::test_that("Phase 0 input validation reports precise failures", {
  X <- data.frame(a = c(0, 1), b = c(1, 2))
  M <- data.frame(a = c(2, 2), b = c(2, 3))
  info <- validate_binomial_data(X, M)
  testthat::expect_identical(info$n, 2L)
  testthat::expect_identical(info$d, 2L)
  testthat::expect_identical(info$coordinate_names, c("a", "b"))

  bad <- as.matrix(X)
  bad[2, 1] <- 1.5
  testthat::expect_error(validate_binomial_data(bad, as.matrix(M)),
                         "row 2, column 1.*value=1.5.*integer-valued")
  bad <- as.matrix(M)
  bad[1, 2] <- 1
  testthat::expect_error(validate_binomial_data(as.matrix(X), bad),
                         "row 1, column 2.*value=1.*at least 2")
  bad <- as.matrix(X)
  bad[2, 2] <- 4
  testthat::expect_error(validate_binomial_data(bad, as.matrix(M)),
                         "row 2, column 2.*value=4.*must not exceed")
  bad_df <- data.frame(a = c(1, 0), b = c("1", "2"))
  testthat::expect_error(validate_binomial_data(bad_df, M), "column 2.*numeric")
})

testthat::test_that("Phase 0 class interfaces and exact empirical centering work", {
  X <- cbind(first = c(1, 3, 2, 4, 1), second = c(2, 1, 4, 3, 2))
  M <- cbind(first = c(5, 6, 5, 7, 4), second = c(6, 5, 7, 6, 5))
  fit <- fit_latent_binomial_cov(X, M)
  testthat::expect_s3_class(fit, "latent_binomial_cov")
  testthat::expect_lt(max(abs(apply(fit$Psi_hat, c(2, 3), mean))), 1e-12)
  testthat::expect_equal(coef(fit, "covariance"), fit$Sigma_stable)
  testthat::expect_equal(coef(fit, "correlation"), fit$R_latent)
  testthat::expect_equal(coef(fit, "partial"), fit$Partial)
  testthat::expect_false("Psi_hat" %in% names(summary(fit)))
  no_if <- fit_latent_binomial_cov(X, M, keep_influence = FALSE)
  testthat::expect_false("Psi_hat" %in% names(no_if))
})
