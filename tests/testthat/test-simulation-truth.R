run_unbiasedness_check <- function(scenario_id, stream, B = 10000L, n = 500L) {
  estimates <- with_rng_stream(stream, {
    out <- matrix(NA_real_, B, 3L, dimnames = list(NULL, c("Sigma11", "Sigma12", "Sigma22")))
    for (b in seq_len(B)) {
      dat <- generate_pairwise_data(n, scenario_id)
      raw <- raw_covariance_estimate(dat$X, dat$M)
      out[b, ] <- c(raw[1, 1], raw[1, 2], raw[2, 2])
    }
    out
  })
  truth_matrix <- pairwise_truth(scenario_id)$Sigma
  truth <- c(Sigma11 = truth_matrix[1, 1], Sigma12 = truth_matrix[1, 2], Sigma22 = truth_matrix[2, 2])
  average <- colMeans(estimates)
  mcse <- apply(estimates, 2L, sd) / sqrt(B)
  difference <- average - truth
  z_mc <- difference / mcse
  list(scenario_id = scenario_id, B = B, n = n, average = average, truth = truth,
       difference = difference, mcse = mcse, z_mcse = z_mc,
       pass = all(abs(difference) <= 4 * mcse))
}

testthat::test_that("T4 raw covariance is unbiased under non-informative denominators", {
  diagnostic <- run_unbiasedness_check("S1", core_test_streams[[4L]])
  testthat::expect_true(diagnostic$pass,
                        info = paste("MCSE ratios:", paste(format(diagnostic$z_mcse), collapse = ", ")))
  cat("T4 MCSE ratios:", paste(sprintf("%s=%.6f", names(diagnostic$z_mcse), diagnostic$z_mcse), collapse = ", "), "\n")
  record_phase1_diagnostic("T4", diagnostic)
})

testthat::test_that("T5 raw covariance is unbiased with informative correlated denominators", {
  diagnostic <- run_unbiasedness_check("S3", core_test_streams[[5L]])
  testthat::expect_true(diagnostic$pass,
                        info = paste("MCSE ratios:", paste(format(diagnostic$z_mcse), collapse = ", ")))
  cat("T5 MCSE ratios:", paste(sprintf("%s=%.6f", names(diagnostic$z_mcse), diagnostic$z_mcse), collapse = ", "), "\n")
  record_phase1_diagnostic("T5", diagnostic)
})

testthat::test_that("T6 simulation generators reproduce analytic latent truth", {
  pairwise_draws <- 20000000L
  matrix_draws <- 2000000L
  diagnostic <- with_rng_stream(core_test_streams[[6L]], {
    scenarios <- read_pairwise_scenarios()
    pairwise <- data.frame(scenario_id = scenarios$scenario_id,
                           empirical = NA_real_, truth = NA_real_, error = NA_real_)
    for (i in seq_len(nrow(scenarios))) {
      s <- scenarios[i, , drop = FALSE]
      U <- rbeta(pairwise_draws, s$a, s$b)
      V <- rbeta(pairwise_draws, s$a, s$b)
      P2 <- if (s$latent_sign == "positive") {
        (1 - s$delta) * U + s$delta * V
      } else {
        (1 - s$delta) * (1 - U) + s$delta * V
      }
      empirical <- cor(U, P2)
      truth <- pairwise_truth(s)$rho
      pairwise$empirical[i] <- empirical
      pairwise$truth[i] <- truth
      pairwise$error[i] <- empirical - truth
      rm(U, V, P2)
      invisible(gc(FALSE))
    }

    latent <- generate_matrix_latent(matrix_draws)
    P <- latent$P
    truth_m <- matrix_truth()
    empirical_mean <- colMeans(P)
    empirical_cov <- cov(P)
    mean_mcse <- apply(P, 2L, sd) / sqrt(matrix_draws)
    covariance_mcse <- matrix(NA_real_, 5L, 5L)
    for (j in seq_len(5L)) {
      for (k in seq_len(5L)) {
        contribution <- (P[, j] - truth_m$mu[j]) * (P[, k] - truth_m$mu[k])
        covariance_mcse[j, k] <- sd(contribution) / sqrt(matrix_draws)
      }
    }
    list(pairwise = pairwise, pairwise_draws = pairwise_draws,
         matrix_draws = matrix_draws,
         matrix_empirical_mean = empirical_mean, matrix_truth_mean = truth_m$mu,
         matrix_mean_mcse = mean_mcse,
         matrix_empirical_cov = empirical_cov, matrix_truth_cov = truth_m$Sigma,
         matrix_cov_mcse = covariance_mcse,
         matrix_mean_z = (empirical_mean - truth_m$mu) / mean_mcse,
         matrix_cov_z = (empirical_cov - truth_m$Sigma) / covariance_mcse)
  })
  testthat::expect_lte(max(abs(diagnostic$pairwise$error)), 5e-4)
  testthat::expect_true(all(abs(diagnostic$matrix_mean_z) <= 4))
  testthat::expect_true(all(abs(diagnostic$matrix_cov_z) <= 4))
  cat(sprintf("T6 max pairwise correlation error: %.17g; max matrix mean MCSE ratio: %.6f; max matrix covariance MCSE ratio: %.6f\n",
              max(abs(diagnostic$pairwise$error)),
              max(abs(diagnostic$matrix_mean_z)),
              max(abs(diagnostic$matrix_cov_z))))
  record_phase1_diagnostic("T6", diagnostic)
})
