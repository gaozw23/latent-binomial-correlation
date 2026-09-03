start_all <- Sys.time()
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

workers <- as.integer(Sys.getenv("PILOT_WORKERS", unset = "4"))
workers <- max(1L, min(workers, parallel::detectCores(logical = TRUE)))
stream_dir <- "config/rng_streams"
streams <- readRDS(file.path(stream_dir, "pairwise_pilot_streams.rds"))
mapping <- read.csv(file.path(stream_dir, "pairwise_pilot_stream_map.csv"), stringsAsFactors = FALSE)
tasks <- canonical_pairwise_tasks(streams, mapping)

dir.create("results/pilot/pairwise", recursive = TRUE, showWarnings = FALSE)
dir.create("results/pilot/checkpoints", recursive = TRUE, showWarnings = FALSE)
cluster <- open_pilot_cluster(workers)
on.exit(close_pilot_cluster(cluster), add = TRUE)

runtime <- list()
shards <- list()
cells <- unique(mapping[c("scenario_id", "n")])
cells <- cells[order(cells$scenario_id, cells$n), , drop = FALSE]
for (cell_id in seq_len(nrow(cells))) {
  scenario_id <- cells$scenario_id[cell_id]
  n <- cells$n[cell_id]
  cell_tasks <- tasks[vapply(tasks, function(task) task$scenario_id == scenario_id && task$n == n, logical(1))]
  final_path <- file.path("results/pilot/pairwise", sprintf("%s_n%d.rds", scenario_id, n))
  checkpoint_path <- file.path("results/pilot/checkpoints", sprintf("%s_n%d_checkpoint.rds", scenario_id, n))
  cell_start <- Sys.time()
  if (file.exists(final_path)) {
    result <- list(results = readRDS(final_path), complete = TRUE,
                   completed = length(cell_tasks), planned = length(cell_tasks))
    resumed <- TRUE
  } else {
    resumed <- file.exists(checkpoint_path)
    result <- run_checkpointed_tasks(
      cell_tasks, "pairwise", checkpoint_path, final_path, cluster,
      checkpoint_every = 25L, max_new_tasks = Inf, load_balanced = TRUE
    )
  }
  if (!result$complete) stop("Pairwise pilot cell did not complete: ", scenario_id, " n=", n, call. = FALSE)
  validate_pairwise_pilot_output(result$results, PAIRWISE_PILOT_REPLICATES)
  if (file.exists(checkpoint_path)) file.remove(checkpoint_path)
  elapsed <- as.numeric(difftime(Sys.time(), cell_start, units = "secs"))
  runtime[[cell_id]] <- data.frame(
    scenario_id = scenario_id, n = n, planned = length(cell_tasks),
    completed = nrow(result$results), failures = sum(!result$results$success),
    warnings = sum(result$results$warning_count),
    nonfinite = sum(result$results$nonfinite_count),
    negative_variance_rate = mean(result$results$negative_raw_variance),
    projection_rate = mean(result$results$projection_active),
    floor_rate = mean(result$results$floor_active),
    elapsed_seconds = elapsed, result_bytes = as.numeric(object.size(result$results)),
    workers = workers, resumed_existing = resumed, stringsAsFactors = FALSE
  )
  shards[[cell_id]] <- result$results
  cat(sprintf("Completed %s n=%d: %d/%d in %.3f seconds\n",
              scenario_id, n, nrow(result$results), length(cell_tasks), elapsed))
}

all_results <- canonical_pairwise_results(do.call(rbind, shards))
validate_pairwise_pilot_output(all_results, 4800L)
atomic_save_rds(all_results, "results/pilot/pairwise_pilot_all.rds")
write.csv(do.call(rbind, runtime), "results/pilot/pairwise_runtime.csv", row.names = FALSE, quote = TRUE)

rmse_streams <- readRDS(file.path(stream_dir, "rmse_pilot_streams.rds"))
rmse_mapping <- read.csv(file.path(stream_dir, "rmse_pilot_stream_map.csv"), stringsAsFactors = FALSE)
summary <- summarise_pairwise_pilot(all_results, rmse_streams, rmse_mapping)
write.csv(summary, "results/pilot/pairwise_summary.csv", row.names = FALSE, quote = TRUE)
atomic_save_rds(summary, "results/pilot/pairwise_summary.rds")

memory <- as.data.frame(gc())
memory$pool <- rownames(memory)
write.csv(memory, "results/pilot/pairwise_memory.csv", row.names = FALSE, quote = TRUE)
cat("Pairwise pilot completed in",
    as.numeric(difftime(Sys.time(), start_all, units = "secs")), "seconds using", workers, "workers\n")
