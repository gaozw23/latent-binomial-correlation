start_all <- Sys.time()
start_utc <- format(start_all, tz = "UTC", usetz = TRUE)
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

workers <- as.integer(Sys.getenv("PHASE3_PAIRWISE_WORKERS", unset = "4"))
workers <- max(1L, min(workers, parallel::detectCores(logical = TRUE)))
stream_dir <- "config/rng_streams"
streams <- readRDS(file.path(stream_dir, "pairwise_full_streams.rds"))
mapping <- read.csv(file.path(stream_dir, "pairwise_full_stream_map.csv"), stringsAsFactors = FALSE)
tasks <- canonical_pairwise_tasks(streams, mapping)
if (length(tasks) != 120000L) stop("Full pairwise task mapping must contain 120000 tasks", call. = FALSE)

dir.create("results/raw/pairwise", recursive = TRUE, showWarnings = FALSE)
dir.create("results/raw/checkpoints/pairwise", recursive = TRUE, showWarnings = FALSE)
dir.create("results/summaries", recursive = TRUE, showWarnings = FALSE)
cluster <- open_phase3_cluster(workers)
on.exit(close_pilot_cluster(cluster), add = TRUE)

cells <- unique(mapping[c("scenario_id", "n")])
cells <- cells[order(cells$scenario_id, cells$n), , drop = FALSE]
runtime <- vector("list", nrow(cells))
shards <- vector("list", nrow(cells))
for (cell_id in seq_len(nrow(cells))) {
  scenario_id <- cells$scenario_id[cell_id]
  n <- cells$n[cell_id]
  cell_tasks <- tasks[vapply(
    tasks, function(task) task$scenario_id == scenario_id && task$n == n, logical(1)
  )]
  final_path <- file.path("results/raw/pairwise", sprintf("%s_n%d.rds", scenario_id, n))
  checkpoint_path <- file.path(
    "results/raw/checkpoints/pairwise", sprintf("%s_n%d_checkpoint.rds", scenario_id, n)
  )
  resumed <- file.exists(checkpoint_path)
  cell_start <- Sys.time()
  result <- run_checkpointed_tasks(
    cell_tasks, "pairwise", checkpoint_path, final_path, cluster,
    checkpoint_every = 500L, max_new_tasks = Inf, load_balanced = TRUE
  )
  if (!result$complete) stop("Full pairwise cell did not complete: ", scenario_id, " n=", n, call. = FALSE)
  validate_pairwise_pilot_output(result$results, PAIRWISE_FULL_REPLICATES)
  if (!identical(result$results, readRDS(checkpoint_path)) ||
      sha256_file(checkpoint_path) != sha256_file(final_path)) {
    stop("Full pairwise checkpoint/final mismatch: ", scenario_id, " n=", n, call. = FALSE)
  }
  elapsed <- as.numeric(difftime(Sys.time(), cell_start, units = "secs"))
  runtime[[cell_id]] <- data.frame(
    scenario_id = scenario_id, n = n, planned = length(cell_tasks),
    completed = nrow(result$results), failures = sum(!result$results$success),
    warnings = sum(result$results$warning_count),
    unexplained_nonfinite = sum(result$results$nonfinite_count),
    expected_undefined_raw = sum(!result$results$rho_raw_valid),
    negative_variance_rate = mean(result$results$negative_raw_variance),
    projection_rate = mean(result$results$projection_active),
    floor_rate = mean(result$results$floor_active),
    elapsed_seconds = elapsed, throughput_per_second = nrow(result$results) / elapsed,
    result_bytes = as.numeric(object.size(result$results)), workers = workers,
    resumed_existing_checkpoint = resumed,
    checkpoint_sha256 = sha256_file(checkpoint_path),
    final_sha256 = sha256_file(final_path), stringsAsFactors = FALSE
  )
  shards[[cell_id]] <- result$results
  cat(sprintf(
    "Completed %s n=%d: %d/%d in %.3f seconds (%.1f replicates/s)\n",
    scenario_id, n, nrow(result$results), length(cell_tasks), elapsed,
    nrow(result$results) / elapsed
  ))
}

all_results <- canonical_pairwise_results(do.call(rbind, shards))
validate_pairwise_pilot_output(all_results, 120000L)
atomic_save_rds(all_results, "results/raw/pairwise_full_all.rds")
runtime <- do.call(rbind, runtime)
write.csv(runtime, "results/summaries/pairwise_full_runtime.csv", row.names = FALSE, quote = TRUE)

rmse_streams <- readRDS(file.path(stream_dir, "rmse_full_streams.rds"))
rmse_mapping <- read.csv(file.path(stream_dir, "rmse_full_stream_map.csv"), stringsAsFactors = FALSE)
summary <- summarise_pairwise_full(all_results, rmse_streams, rmse_mapping)
atomic_save_rds(summary, "results/summaries/pairwise_full_summary.rds")
write.csv(summary, "results/summaries/pairwise_full_summary.csv", row.names = FALSE, quote = TRUE)
diagnostics <- summarise_pairwise_diagnostics(all_results)
atomic_save_rds(diagnostics, "results/summaries/pairwise_full_diagnostics.rds")
write.csv(diagnostics, "results/summaries/pairwise_full_diagnostics.csv", row.names = FALSE, quote = TRUE)

memory <- as.data.frame(gc())
memory$pool <- rownames(memory)
write.csv(memory, "results/summaries/pairwise_full_memory.csv", row.names = FALSE, quote = TRUE)
end_all <- Sys.time()
metadata <- list(
  start_utc = start_utc, end_utc = format(end_all, tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(end_all, start_all, units = "secs")),
  workers = workers, tasks = nrow(all_results),
  throughput_per_second = nrow(all_results) /
    as.numeric(difftime(end_all, start_all, units = "secs")),
  approximate_memory = memory
)
atomic_save_rds(metadata, "results/summaries/pairwise_full_execution_metadata.rds")
cat("Full pairwise simulation completed in", metadata$elapsed_seconds,
    "seconds using", workers, "workers\n")
