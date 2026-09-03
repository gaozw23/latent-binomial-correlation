start_all <- Sys.time()
start_utc <- format(start_all, tz = "UTC", usetz = TRUE)
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

workers <- as.integer(Sys.getenv("PHASE4_MATRIX_WORKERS", unset = "4"))
workers <- max(1L, min(workers, parallel::detectCores(logical = TRUE)))
streams <- readRDS("config/rng_streams/matrix_full_streams.rds")
mapping <- read.csv("config/rng_streams/matrix_full_stream_map.csv", stringsAsFactors = FALSE)
tasks <- canonical_matrix_full_tasks(streams, mapping)
if (length(tasks) != 15000L) stop("Full matrix mapping must contain 15000 tasks", call. = FALSE)

dir.create("results/raw/matrix", recursive = TRUE, showWarnings = FALSE)
dir.create("results/raw/checkpoints/matrix", recursive = TRUE, showWarnings = FALSE)
dir.create("results/summaries", recursive = TRUE, showWarnings = FALSE)
cluster <- open_phase3_cluster(workers)
on.exit(close_pilot_cluster(cluster), add = TRUE)

runtime <- vector("list", length(MATRIX_N_VALUES))
shards <- vector("list", length(MATRIX_N_VALUES))
for (cell_id in seq_along(MATRIX_N_VALUES)) {
  n <- MATRIX_N_VALUES[cell_id]
  cell_tasks <- tasks[vapply(tasks, function(task) task$n == n, logical(1))]
  final_path <- file.path("results/raw/matrix", sprintf("matrix_n%d.rds", n))
  checkpoint_path <- file.path("results/raw/checkpoints/matrix", sprintf("matrix_n%d_checkpoint.rds", n))
  resumed <- file.exists(checkpoint_path)
  cell_start <- Sys.time()
  result <- run_checkpointed_matrix_full_tasks(
    cell_tasks, checkpoint_path, final_path, cluster,
    checkpoint_every = 500L, load_balanced = TRUE
  )
  if (!result$complete) stop("Full matrix cell did not complete: n=", n, call. = FALSE)
  validate_matrix_full_output(result$results, 5000L)
  if (!identical(result$results, readRDS(checkpoint_path)) ||
      sha256_file(checkpoint_path) != sha256_file(final_path)) {
    stop("Full matrix checkpoint/final mismatch: n=", n, call. = FALSE)
  }
  elapsed <- as.numeric(difftime(Sys.time(), cell_start, units = "secs"))
  runtime[[cell_id]] <- data.frame(
    n = n, planned = 5000L, completed = nrow(result$results),
    failures = sum(!result$results$success), warnings = sum(result$results$warning_count),
    unexplained_nonfinite = sum(result$results$nonfinite_count),
    theorem_violations = sum(!result$results$projection_inequality_pass),
    elapsed_seconds = elapsed, throughput_per_second = nrow(result$results) / elapsed,
    result_bytes = as.numeric(object.size(result$results)), workers = workers,
    resumed_existing_checkpoint = resumed,
    checkpoint_sha256 = sha256_file(checkpoint_path), final_sha256 = sha256_file(final_path),
    stringsAsFactors = FALSE
  )
  shards[[cell_id]] <- result$results
  cat(sprintf("Completed matrix n=%d: 5000/5000 in %.3f seconds (%.1f replicates/s)\n",
              n, elapsed, nrow(result$results) / elapsed))
}

all_results <- canonical_matrix_full_results(do.call(rbind, shards))
validate_matrix_full_output(all_results, 15000L)
atomic_save_rds(all_results, "results/raw/matrix_full_all.rds")
runtime <- do.call(rbind, runtime)
write.csv(runtime, "results/summaries/matrix_full_runtime.csv", row.names = FALSE, quote = TRUE)
summary <- summarise_matrix_full(all_results)
atomic_save_rds(summary, "results/summaries/matrix_full_summary.rds")
write.csv(summary, "results/summaries/matrix_full_summary.csv", row.names = FALSE, quote = TRUE)
correlation_pairs <- summarise_matrix_pairs(all_results, "correlation")
partial_pairs <- summarise_matrix_pairs(all_results, "partial")
atomic_save_rds(correlation_pairs, "results/summaries/matrix_correlation_pair_summary.rds")
write.csv(correlation_pairs, "results/summaries/matrix_correlation_pair_summary.csv",
          row.names = FALSE, quote = TRUE)
atomic_save_rds(partial_pairs, "results/summaries/matrix_partial_pair_summary.rds")
write.csv(partial_pairs, "results/summaries/matrix_partial_pair_summary.csv",
          row.names = FALSE, quote = TRUE)
truth <- matrix_truth()
atomic_save_rds(truth, "results/summaries/matrix_analytic_truth.rds")

memory <- as.data.frame(gc())
memory$pool <- rownames(memory)
write.csv(memory, "results/summaries/matrix_full_memory.csv", row.names = FALSE, quote = TRUE)
end_all <- Sys.time()
metadata <- list(
  start_utc = start_utc, end_utc = format(end_all, tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(end_all, start_all, units = "secs")),
  workers = workers, tasks = nrow(all_results),
  throughput_per_second = nrow(all_results) /
    as.numeric(difftime(end_all, start_all, units = "secs")),
  approximate_memory = memory
)
atomic_save_rds(metadata, "results/summaries/matrix_full_execution_metadata.rds")
cat("Full matrix simulation completed in", metadata$elapsed_seconds,
    "seconds using", workers, "workers\n")
