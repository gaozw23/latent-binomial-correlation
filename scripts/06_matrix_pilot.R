start_all <- Sys.time()
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

workers <- as.integer(Sys.getenv("PILOT_WORKERS", unset = "4"))
workers <- max(1L, min(workers, parallel::detectCores(logical = TRUE)))
stream_dir <- "config/rng_streams"
streams <- readRDS(file.path(stream_dir, "matrix_pilot_streams.rds"))
mapping <- read.csv(file.path(stream_dir, "matrix_pilot_stream_map.csv"), stringsAsFactors = FALSE)
tasks <- canonical_matrix_tasks(streams, mapping)

dir.create("results/pilot/matrix", recursive = TRUE, showWarnings = FALSE)
dir.create("results/pilot/checkpoints", recursive = TRUE, showWarnings = FALSE)
cluster <- open_pilot_cluster(workers)
on.exit(close_pilot_cluster(cluster), add = TRUE)

runtime <- list()
shards <- list()
for (cell_id in seq_along(MATRIX_N_VALUES)) {
  n <- MATRIX_N_VALUES[cell_id]
  cell_tasks <- tasks[vapply(tasks, function(task) task$n == n, logical(1))]
  final_path <- file.path("results/pilot/matrix", sprintf("matrix_n%d.rds", n))
  checkpoint_path <- file.path("results/pilot/checkpoints", sprintf("matrix_n%d_checkpoint.rds", n))
  cell_start <- Sys.time()
  if (file.exists(final_path)) {
    result <- list(results = readRDS(final_path), complete = TRUE,
                   completed = length(cell_tasks), planned = length(cell_tasks))
    resumed <- TRUE
  } else {
    resumed <- file.exists(checkpoint_path)
    result <- run_checkpointed_tasks(
      cell_tasks, "matrix", checkpoint_path, final_path, cluster,
      checkpoint_every = 25L, max_new_tasks = Inf, load_balanced = TRUE
    )
  }
  if (!result$complete) stop("Matrix pilot cell did not complete: n=", n, call. = FALSE)
  validate_matrix_pilot_output(result$results, MATRIX_PILOT_REPLICATES)
  if (file.exists(checkpoint_path)) file.remove(checkpoint_path)
  elapsed <- as.numeric(difftime(Sys.time(), cell_start, units = "secs"))
  runtime[[cell_id]] <- data.frame(
    n = n, planned = length(cell_tasks), completed = nrow(result$results),
    failures = sum(!result$results$success), warnings = sum(result$results$warning_count),
    nonfinite = sum(result$results$nonfinite_count),
    theorem_violations = sum(!result$results$projection_inequality_pass),
    minimum_theorem_margin = min(result$results$projection_inequality_margin),
    projection_rate = mean(result$results$projection_active),
    floor_rate = mean(result$results$floor_active),
    elapsed_seconds = elapsed, result_bytes = as.numeric(object.size(result$results)),
    workers = workers, resumed_existing = resumed, stringsAsFactors = FALSE
  )
  shards[[cell_id]] <- result$results
  cat(sprintf("Completed matrix n=%d: %d/%d in %.3f seconds\n",
              n, nrow(result$results), length(cell_tasks), elapsed))
}

all_results <- canonical_matrix_results(do.call(rbind, shards))
validate_matrix_pilot_output(all_results, 1500L)
atomic_save_rds(all_results, "results/pilot/matrix_pilot_all.rds")
write.csv(do.call(rbind, runtime), "results/pilot/matrix_runtime.csv", row.names = FALSE, quote = TRUE)
summary <- summarise_matrix_pilot(all_results)
write.csv(summary, "results/pilot/matrix_summary.csv", row.names = FALSE, quote = TRUE)
atomic_save_rds(summary, "results/pilot/matrix_summary.rds")
memory <- as.data.frame(gc())
memory$pool <- rownames(memory)
write.csv(memory, "results/pilot/matrix_memory.csv", row.names = FALSE, quote = TRUE)
cat("Matrix pilot completed in",
    as.numeric(difftime(Sys.time(), start_all, units = "secs")), "seconds using", workers, "workers\n")
