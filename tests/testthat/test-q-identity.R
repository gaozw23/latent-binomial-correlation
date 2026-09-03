testthat::test_that("T1 exact enumeration verifies the Q identity", {
  max_error <- 0
  worst <- NULL
  for (m in 2:20) {
    x <- 0:m
    q <- x * (m - x) / (m^2 * (m - 1))
    for (p in seq(0, 1, by = 0.025)) {
      lhs <- sum(dbinom(x, m, p) * q)
      rhs <- p * (1 - p) / m
      error <- abs(lhs - rhs)
      if (error > max_error) {
        max_error <- error
        worst <- c(m = m, p = p, lhs = lhs, rhs = rhs)
      }
      testthat::expect_lte(error, 5e-15)
    }
  }
  cat(sprintf("T1 max absolute error: %.17g\n", max_error))
  record_phase1_diagnostic("T1", list(max_absolute_error = max_error, worst_case = worst,
                                       tolerance = 5e-15, enumerated_cases = 19L * 41L))
})
