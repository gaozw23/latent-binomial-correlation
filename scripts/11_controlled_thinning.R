start_time <- Sys.time()
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

dir.create("results/summaries", recursive = TRUE, showWarnings = FALSE)
dir.create("results/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("config/rng_streams", recursive = TRUE, showWarnings = FALSE)

data <- read_frozen_aec()
caps <- c(5L, 10L, 20L, 50L, 100L)
streams <- make_streams(length(caps) * 2000L, MASTER_SEEDS[["aec_thinning"]])
saveRDS(streams, "config/rng_streams/aec_thinning_streams.rds", compress = "xz")

result <- controlled_binomial_thinning(data, streams, caps, 2000L)
write.csv(result$mapping, "config/rng_streams/aec_thinning_stream_map.csv",
          row.names = FALSE, quote = TRUE)
saveRDS(result$replicates, "results/raw/aec_thinning.rds", compress = "xz")
write.csv(result$summary, "results/summaries/aec_thinning_summary.csv",
          row.names = FALSE, quote = TRUE)
saveRDS(list(
  started_utc = format(start_time, tz = "UTC", usetz = TRUE),
  ended_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
  tasks = nrow(result$replicates), caps = caps
), "results/summaries/aec_thinning_execution_metadata.rds", compress = "xz")
cat(sprintf("Controlled thinning complete: %d tasks; failures/non-finite=%d\n",
            nrow(result$replicates), sum(!result$replicates$valid)))
