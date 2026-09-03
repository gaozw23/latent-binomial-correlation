start_time <- Sys.time()
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
dir.create("results/audits", recursive = TRUE, showWarnings = FALSE)

source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

log_file <- "results/audits/phase1_test_output.log"
zz <- file(log_file, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({
  while (sink.number(type = "message") > 2L) sink(type = "message")
  while (sink.number(type = "output") > 0L) sink(type = "output")
  close(zz)
}, add = TRUE)

cat("Phase 1 tests started UTC:", format(start_time, tz = "UTC", usetz = TRUE), "\n")
results <- testthat::test_dir(
  "tests/testthat",
  reporter = "summary",
  stop_on_failure = TRUE,
  stop_on_warning = FALSE
)
saveRDS(results, "results/audits/phase1_testthat_results.rds")
cat("Phase 1 tests ended UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
cat("Elapsed seconds:", as.numeric(difftime(Sys.time(), start_time, units = "secs")), "\n")
