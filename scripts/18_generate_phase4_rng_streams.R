start <- Sys.time()
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

manifest <- generate_phase4_stream_files()
dir.create("results/manifests", recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, "results/manifests/phase4_rng_stream_manifest.csv", row.names = FALSE, quote = TRUE)
print(manifest)
cat("Matrix-simulation RNG stream generation completed in",
    as.numeric(difftime(Sys.time(), start, units = "secs")), "seconds\n")
