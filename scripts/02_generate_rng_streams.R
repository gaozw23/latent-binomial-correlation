start <- Sys.time()
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)
manifest <- generate_phase2_stream_files()
manifest$master_seed <- c(MASTER_SEEDS[["pairwise"]], NA, MASTER_SEEDS[["matrix"]], NA,
                          MASTER_SEEDS[["rmse_mcse"]], NA)
write.csv(manifest, "results/pilot/rng_stream_manifest.csv", row.names = FALSE, quote = TRUE)
cat("Canonical Phase 2 RNG streams generated and validated in",
    as.numeric(difftime(Sys.time(), start, units = "secs")), "seconds\n")
print(manifest)
