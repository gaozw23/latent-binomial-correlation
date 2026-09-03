testthat::test_that("T3 PSD projection has the required geometry and loss inequality", {
  diagnostics <- with_rng_stream(core_test_streams[[3L]], {
    max_symmetry_error <- 0
    min_projected_eigenvalue <- Inf
    max_distance_identity_error <- 0
    min_loss_margin <- Inf
    input_count <- 0L
    target_count <- 0L
    for (d in c(2L, 5L, 10L)) {
      for (replicate in seq_len(100L)) {
        Z <- matrix(rnorm(d * d), d, d)
        A <- symmetrise_matrix(Z)
        projection <- project_psd(A)
        input_count <- input_count + 1L
        symmetry_error <- norm(projection$matrix - t(projection$matrix), "F")
        minimum <- min(eigen(projection$matrix, symmetric = TRUE, only.values = TRUE)$values)
        expected_distance <- sqrt(sum(pmin(projection$raw_eigenvalues, 0)^2))
        distance_error <- abs(projection$frobenius_adjustment - expected_distance)
        testthat::expect_lte(symmetry_error, 1e-14)
        testthat::expect_gte(minimum, -1e-12)
        testthat::expect_lte(distance_error, 1e-12 * (1 + expected_distance))
        max_symmetry_error <- max(max_symmetry_error, symmetry_error)
        min_projected_eigenvalue <- min(min_projected_eigenvalue, minimum)
        max_distance_identity_error <- max(max_distance_identity_error, distance_error)
        for (target_id in seq_len(20L)) {
          B <- matrix(rnorm(d * d), d, d)
          target <- tcrossprod(B) / d
          raw_loss <- norm(A - target, "F")
          projected_loss <- norm(projection$matrix - target, "F")
          margin <- raw_loss + 1e-10 * (1 + raw_loss) - projected_loss
          testthat::expect_gte(margin, 0)
          min_loss_margin <- min(min_loss_margin, margin)
          target_count <- target_count + 1L
        }
      }
    }
    list(input_count = input_count, target_count = target_count,
         max_symmetry_error = max_symmetry_error,
         min_projected_eigenvalue = min_projected_eigenvalue,
         max_distance_identity_error = max_distance_identity_error,
         min_loss_margin = min_loss_margin)
  })
  cat(sprintf(paste0("T3 inputs: %d; PSD targets: %d; minimum eigenvalue: %.17g; ",
                     "max distance error: %.17g; minimum loss margin: %.17g\n"),
              diagnostics$input_count, diagnostics$target_count,
              diagnostics$min_projected_eigenvalue,
              diagnostics$max_distance_identity_error,
              diagnostics$min_loss_margin))
  record_phase1_diagnostic("T3", diagnostics)
})
