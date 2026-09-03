testthat::test_that("T12 worker-count reproducibility is prepared for Phase 2", {
  if (!identical(Sys.getenv("RUN_PHASE2_T12"), "true")) {
    testthat::skip("T12 requires Phase 2 authorization; set RUN_PHASE2_T12=true")
  }
  stream_dir <- file.path(PROJECT_ROOT, "config", "rng_streams")
  streams <- readRDS(file.path(stream_dir, "pairwise_pilot_streams.rds"))
  mapping <- read.csv(file.path(stream_dir, "pairwise_pilot_stream_map.csv"), stringsAsFactors = FALSE)
  tasks <- canonical_pairwise_tasks(streams, mapping)
  diagnostic <- run_t12_reproducibility(tasks, file.path(PROJECT_ROOT, "results", "audits"))
  testthat::expect_true(diagnostic$pass)
  testthat::expect_identical(diagnostic$workers_1_hash, diagnostic$workers_2_hash)
  record_phase1_diagnostic("T12", diagnostic)
})
