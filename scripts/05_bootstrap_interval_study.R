start_all <- Sys.time()
start_utc <- format(start_all, tz = "UTC", usetz = TRUE)
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

workers <- as.integer(Sys.getenv("PHASE3_BOOTSTRAP_WORKERS", unset = "4"))
workers <- max(1L, min(workers, parallel::detectCores(logical = TRUE)))

sanity_start <- Sys.time()
sanity_cluster <- open_phase3_cluster(min(2L, workers))
sanity_ok <- if (is.null(sanity_cluster)) TRUE else all(unlist(parallel::clusterEvalQ(sanity_cluster, 1L + 1L == 2L)))
close_pilot_cluster(sanity_cluster)
if (!sanity_ok) stop("Bootstrap PSOCK worker startup sanity check failed", call. = FALSE)
sanity_seconds <- as.numeric(difftime(Sys.time(), sanity_start, units = "secs"))
cat(sprintf("Bootstrap worker startup sanity check passed in %.3f seconds\n", sanity_seconds))

stream_dir <- "config/rng_streams"
streams <- readRDS(file.path(stream_dir, "interval_bootstrap_streams.rds"))
mapping <- read.csv(file.path(stream_dir, "interval_bootstrap_stream_map.csv"), stringsAsFactors = FALSE)
tasks <- canonical_bootstrap_tasks(streams, mapping)
if (length(tasks) != 4000L) stop("Focused bootstrap mapping must contain 4000 outer tasks", call. = FALSE)

dir.create("results/raw/bootstrap", recursive = TRUE, showWarnings = FALSE)
dir.create("results/raw/checkpoints/bootstrap", recursive = TRUE, showWarnings = FALSE)
dir.create("results/summaries", recursive = TRUE, showWarnings = FALSE)
cluster <- open_phase3_cluster(workers)
on.exit(close_pilot_cluster(cluster), add = TRUE)

cells <- unique(mapping[c("scenario_id", "n")])
cells <- cells[order(cells$scenario_id, cells$n), , drop = FALSE]
runtime <- vector("list", nrow(cells))
shards <- vector("list", nrow(cells))
runtime_progress_path <- "results/summaries/bootstrap_focused_runtime_progress.csv"
runtime_progress <- if (file.exists(runtime_progress_path)) {
  read.csv(runtime_progress_path, stringsAsFactors = FALSE)
} else data.frame()
for (cell_id in seq_len(nrow(cells))) {
  scenario_id <- cells$scenario_id[cell_id]
  n <- cells$n[cell_id]
  cell_tasks <- tasks[vapply(
    tasks, function(task) task$scenario_id == scenario_id && task$n == n, logical(1)
  )]
  final_path <- file.path("results/raw/bootstrap", sprintf("%s_n%d.rds", scenario_id, n))
  checkpoint_path <- file.path(
    "results/raw/checkpoints/bootstrap", sprintf("%s_n%d_checkpoint.rds", scenario_id, n)
  )
  prior_runtime <- if (nrow(runtime_progress)) {
    runtime_progress[
      runtime_progress$scenario_id == scenario_id & runtime_progress$n == n, , drop = FALSE
    ]
  } else data.frame()
  if (file.exists(final_path) && file.exists(checkpoint_path) && nrow(prior_runtime) == 1L) {
    existing_result <- readRDS(final_path)
    validate_bootstrap_output(existing_result, 1000L)
    if (!identical(existing_result, readRDS(checkpoint_path)) ||
        sha256_file(checkpoint_path) != sha256_file(final_path)) {
      stop("Existing bootstrap checkpoint/final mismatch: ", scenario_id, " n=", n, call. = FALSE)
    }
    runtime[[cell_id]] <- prior_runtime
    shards[[cell_id]] <- existing_result
    cat(sprintf("Reused verified bootstrap %s n=%d: 1000/1000 outer datasets\n", scenario_id, n))
    next
  }
  resumed <- file.exists(checkpoint_path)
  cell_start <- Sys.time()
  result <- run_checkpointed_bootstrap_tasks(
    cell_tasks, checkpoint_path, final_path, cluster,
    checkpoint_every = 25L, load_balanced = TRUE
  )
  if (!result$complete) stop("Focused bootstrap cell did not complete: ", scenario_id, " n=", n, call. = FALSE)
  validate_bootstrap_output(result$results, 1000L)
  if (!identical(result$results, readRDS(checkpoint_path)) ||
      sha256_file(checkpoint_path) != sha256_file(final_path)) {
    stop("Bootstrap checkpoint/final mismatch: ", scenario_id, " n=", n, call. = FALSE)
  }
  elapsed <- as.numeric(difftime(Sys.time(), cell_start, units = "secs"))
  runtime[[cell_id]] <- data.frame(
    scenario_id = scenario_id, n = n, planned_outer = length(cell_tasks),
    completed_outer = nrow(result$results), failed_outer = sum(!result$results$success),
    inner_draws_attempted = sum(result$results$bootstrap_draws_attempted),
    inner_draws_valid = sum(result$results$bootstrap_draws_valid),
    inner_failures = sum(result$results$bootstrap_failure_count),
    elapsed_seconds = elapsed, outer_tasks_per_second = nrow(result$results) / elapsed,
    inner_draws_per_second = sum(result$results$bootstrap_draws_attempted) / elapsed,
    result_bytes = as.numeric(object.size(result$results)), workers = workers,
    resumed_existing_checkpoint = resumed,
    checkpoint_sha256 = sha256_file(checkpoint_path),
    final_sha256 = sha256_file(final_path), stringsAsFactors = FALSE
  )
  shards[[cell_id]] <- result$results
  completed_runtime <- do.call(rbind, runtime[!vapply(runtime, is.null, logical(1))])
  write.csv(completed_runtime, runtime_progress_path, row.names = FALSE, quote = TRUE)
  cat(sprintf(
    "Completed bootstrap %s n=%d: %d/%d outer datasets in %.3f seconds (%.1f inner draws/s)\n",
    scenario_id, n, nrow(result$results), length(cell_tasks), elapsed,
    sum(result$results$bootstrap_draws_attempted) / elapsed
  ))
}

all_results <- canonical_bootstrap_results(do.call(rbind, shards))
validate_bootstrap_output(all_results, 4000L)
atomic_save_rds(all_results, "results/raw/bootstrap_focused_all.rds")
runtime <- do.call(rbind, runtime)
write.csv(runtime, "results/summaries/bootstrap_focused_runtime.csv", row.names = FALSE, quote = TRUE)
summary <- summarise_focused_bootstrap(all_results)
atomic_save_rds(summary, "results/summaries/bootstrap_focused_summary.rds")
write.csv(summary, "results/summaries/bootstrap_focused_summary.csv", row.names = FALSE, quote = TRUE)

memory <- as.data.frame(gc())
memory$pool <- rownames(memory)
write.csv(memory, "results/summaries/bootstrap_focused_memory.csv", row.names = FALSE, quote = TRUE)
end_all <- Sys.time()
metadata <- list(
  start_utc = start_utc, end_utc = format(end_all, tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(end_all, start_all, units = "secs")),
  workers = workers, outer_tasks = nrow(all_results),
  inner_draws = sum(all_results$bootstrap_draws_attempted),
  outer_tasks_per_second = nrow(all_results) /
    as.numeric(difftime(end_all, start_all, units = "secs")),
  inner_draws_per_second = sum(all_results$bootstrap_draws_attempted) /
    as.numeric(difftime(end_all, start_all, units = "secs")),
  worker_startup_sanity_pass = sanity_ok, worker_startup_sanity_seconds = sanity_seconds,
  approximate_memory = memory
)
atomic_save_rds(metadata, "results/summaries/bootstrap_focused_execution_metadata.rds")
cat("Focused bootstrap completed in", metadata$elapsed_seconds,
    "seconds using", workers, "workers\n")
