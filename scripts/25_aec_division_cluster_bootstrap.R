start_time <- Sys.time()
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

master_seed <- 20260909L
replicates <- 4999L
years <- c(2019L, 2022L)
workers <- as.integer(Sys.getenv("PHASE5C_WORKERS", unset = "4"))
workers <- max(1L, min(workers, parallel::detectCores(logical = TRUE)))

dir.create("config/rng_streams", recursive = TRUE, showWarnings = FALSE)
dir.create("results/raw/checkpoints", recursive = TRUE, showWarnings = FALSE)
dir.create("results/summaries", recursive = TRUE, showWarnings = FALSE)

data <- read_frozen_aec("results/derived/aec_matched_frozen_with_divisions.csv")
if (nrow(data) != 6874L || anyNA(data$DivisionID_2019) || anyNA(data$DivisionID_2022)) {
  stop("The Phase 5C enriched derivative is incomplete", call. = FALSE)
}
structure <- audit_aec_division_structure(data)
if (isTRUE(structure$summary$all_division_ids_identical)) years <- 2019L

mapping <- do.call(rbind, lapply(years, function(year) data.frame(
  stream_id = seq_len(replicates), clustering_year = year,
  replicate = seq_len(replicates), master_seed = master_seed,
  study_component = "aec_state_stratified_division_cluster_bootstrap",
  stringsAsFactors = FALSE
)))
mapping$stream_id <- seq_len(nrow(mapping))
streams <- make_streams(nrow(mapping), master_seed)
saveRDS(streams, "config/rng_streams/aec_division_cluster_bootstrap_streams.rds",
        compress = "xz")
write.csv(mapping, "config/rng_streams/aec_division_cluster_bootstrap_stream_map.csv",
          row.names = FALSE, quote = TRUE)

cluster <- open_phase3_cluster(workers)
on.exit(close_pilot_cluster(cluster), add = TRUE)
summaries <- list()
runtime <- list()
for (year in years) {
  year_map <- mapping[mapping$clustering_year == year, , drop = FALSE]
  tasks <- lapply(seq_len(nrow(year_map)), function(i) list(
    stream_id = as.integer(year_map$stream_id[i]),
    clustering_year = as.integer(year_map$clustering_year[i]),
    replicate = as.integer(year_map$replicate[i]),
    seed = streams[[year_map$stream_id[i]]]
  ))
  cluster_indices <- aec_division_cluster_indices(data, year)
  checkpoint_path <- sprintf(
    "results/raw/checkpoints/aec_division_cluster_bootstrap_%d_checkpoint.rds", year
  )
  final_path <- sprintf("results/raw/aec_division_cluster_bootstrap_%d.rds", year)
  started <- Sys.time()
  result <- run_checkpointed_aec_division_tasks(
    tasks, data, cluster_indices, checkpoint_path, final_path,
    cluster = cluster, checkpoint_every = 250L
  )
  if (!result$complete) stop("Phase 5C bootstrap did not complete for ", year, call. = FALSE)
  if (nrow(result$results) != replicates || anyDuplicated(result$results$replicate) ||
      !identical(result$results$replicate, seq_len(replicates))) {
    stop("Phase 5C bootstrap replicate identities are incomplete for ", year, call. = FALSE)
  }
  summaries[[as.character(year)]] <- summarise_aec_division_bootstrap(
    result$results, data, year, workers, master_seed
  )
  write.csv(
    summaries[[as.character(year)]],
    sprintf("results/summaries/aec_division_cluster_bootstrap_%d.csv", year),
    row.names = FALSE, quote = TRUE
  )
  runtime[[as.character(year)]] <- data.frame(
    clustering_year = year, workers = workers, planned = replicates,
    completed = nrow(result$results), failures = sum(!result$results$success),
    warnings = sum(result$results$warning_count),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    checkpoint_path = checkpoint_path, final_path = final_path,
    stringsAsFactors = FALSE
  )
  cat(sprintf("Completed %d division bootstrap: %d/%d valid=%d\n", year,
              nrow(result$results), replicates, sum(result$results$success)))
}
runtime <- do.call(rbind, runtime)
write.csv(runtime, "results/summaries/aec_division_cluster_bootstrap_runtime.csv",
          row.names = FALSE, quote = TRUE)
saveRDS(list(
  started_utc = format(start_time, tz = "UTC", usetz = TRUE),
  ended_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
  master_seed = master_seed, workers = workers, years = years,
  replicates_per_year = replicates
), "results/summaries/aec_division_cluster_bootstrap_execution_metadata.rds",
compress = "xz")
