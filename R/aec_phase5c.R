.phase5c_required_frozen_columns <- c(
  "match_id", "match_route", "State", "PollingPlaceID_2019",
  "PollingPlaceID_2022", "X_2019", "M_2019", "Y_2019",
  "X_2022", "M_2022", "Y_2022"
)

.phase5c_csv_quote <- function(x) {
  x <- as.character(x)
  paste0('"', gsub('"', '""', x, fixed = TRUE), '"')
}

create_aec_division_enrichment <- function(
    frozen_path = "results/derived/aec_matched_frozen.csv",
    output_path = "results/derived/aec_matched_frozen_with_divisions.csv",
    provenance_path = "results/audits/aec_division_enrichment_provenance.csv") {
  if (file.exists(output_path) || file.exists(provenance_path)) {
    stop("Refusing to overwrite an existing Phase 5C enrichment artifact", call. = FALSE)
  }
  frozen <- read.csv(frozen_path, stringsAsFactors = FALSE, check.names = FALSE,
                     colClasses = "character", na.strings = NULL)
  missing <- setdiff(.phase5c_required_frozen_columns, names(frozen))
  if (length(missing) || nrow(frozen) != 6874L || anyDuplicated(frozen$match_id)) {
    stop("The canonical Phase 5B frozen cohort failed its structural precheck", call. = FALSE)
  }

  inventory <- parse_aec_inventory("results/manifests/aec_acquisition_manifest.csv")
  attached <- list()
  provenance <- list()
  for (year in c(2019L, 2022L)) {
    metadata <- inventory$metadata[[as.character(year)]]
    metadata <- metadata[metadata$polling_place_type_id == "1", , drop = FALSE]
    metadata_key <- paste(metadata$state_ab, metadata$polling_place_id, sep = "|")
    duplicate_key_count <- sum(duplicated(metadata_key))
    if (duplicate_key_count != 0L) {
      stop(year, " official ordinary metadata has duplicate state/PollingPlaceID keys", call. = FALSE)
    }
    id_column <- paste0("PollingPlaceID_", year)
    frozen_key <- paste(frozen$State, frozen[[id_column]], sep = "|")
    index <- match(frozen_key, metadata_key)
    if (anyNA(index) || length(index) != nrow(frozen)) {
      stop(year, " division enrichment has a zero-match frozen row", call. = FALSE)
    }
    if (!identical(frozen$State, metadata$state_ab[index]) ||
        !identical(frozen[[id_column]], metadata$polling_place_id[index])) {
      stop(year, " division enrichment failed exact state/ID consistency", call. = FALSE)
    }
    division_id <- metadata$division_id[index]
    division_name <- metadata$division_name[index]
    if (any(!nzchar(division_id))) stop(year, " DivisionID coverage is incomplete", call. = FALSE)
    attached[[paste0("DivisionID_", year)]] <- division_id
    attached[[paste0("DivisionNm_", year)]] <- division_name

    manifest_row <- inventory$manifest[
      inventory$manifest$election_year == year &
        inventory$manifest$file_role == "polling_places_metadata", , drop = FALSE
    ]
    if (nrow(manifest_row) != 1L) stop("Official metadata provenance is ambiguous", call. = FALSE)
    schema <- inventory$schema_inventory[
      inventory$schema_inventory$election_year == year &
        inventory$schema_inventory$file_role == "polling_places_metadata", , drop = FALSE
    ]
    source_id <- schema$original_column_name[schema$canonical_concept == "division_id"]
    source_name <- schema$original_column_name[schema$canonical_concept == "division_name"]
    if (length(source_id) != 1L || length(source_name) != 1L) {
      stop("Official division source-column provenance is ambiguous", call. = FALSE)
    }
    provenance[[length(provenance) + 1L]] <- data.frame(
      output_field = paste0("DivisionID_", year), election_year = year,
      official_source_file = manifest_row$official_basename,
      official_source_path = manifest_row$local_path,
      source_column_name = source_id,
      join_key = paste0("State + PollingPlaceID_", year),
      row_match_count = nrow(frozen), missing_count = sum(!nzchar(division_id)),
      duplicate_key_count = duplicate_key_count, stringsAsFactors = FALSE
    )
    provenance[[length(provenance) + 1L]] <- data.frame(
      output_field = paste0("DivisionNm_", year), election_year = year,
      official_source_file = manifest_row$official_basename,
      official_source_path = manifest_row$local_path,
      source_column_name = source_name,
      join_key = paste0("State + PollingPlaceID_", year),
      row_match_count = nrow(frozen), missing_count = sum(!nzchar(division_name)),
      duplicate_key_count = duplicate_key_count, stringsAsFactors = FALSE
    )
  }

  # Preserve each original CSV line byte-for-byte as the prefix, appending only new fields.
  original_lines <- readLines(frozen_path, warn = FALSE, encoding = "UTF-8")
  if (length(original_lines) != nrow(frozen) + 1L) {
    stop("Frozen CSV line count is incompatible with safe append-only enrichment", call. = FALSE)
  }
  header_suffix <- paste(c("DivisionID_2019", "DivisionNm_2019",
                           "DivisionID_2022", "DivisionNm_2022"), collapse = ",")
  values <- cbind(
    attached$DivisionID_2019, attached$DivisionNm_2019,
    attached$DivisionID_2022, attached$DivisionNm_2022
  )
  suffix <- apply(values, 1L, function(row) paste(.phase5c_csv_quote(row), collapse = ","))
  enriched_lines <- c(paste0(original_lines[1L], ",", header_suffix),
                      paste0(original_lines[-1L], ",", suffix))
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enriched_lines, output_path, useBytes = TRUE)

  enriched <- read.csv(output_path, stringsAsFactors = FALSE, check.names = FALSE,
                       colClasses = "character", na.strings = NULL)
  if (nrow(enriched) != nrow(frozen) ||
      !identical(enriched[names(frozen)], frozen) ||
      !identical(enriched$match_id, frozen$match_id) || anyDuplicated(enriched$match_id) ||
      any(!nzchar(enriched$DivisionID_2019)) || any(!nzchar(enriched$DivisionID_2022))) {
    stop("Written Phase 5C enrichment failed post-write validation", call. = FALSE)
  }
  derivative_sha256 <- sha256_file(output_path)
  provenance <- do.call(rbind, provenance)
  provenance$enriched_derivative <- output_path
  provenance$enriched_derivative_sha256 <- derivative_sha256
  dir.create(dirname(provenance_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(provenance, provenance_path, row.names = FALSE, quote = TRUE)
  Sys.chmod(output_path, mode = "0444", use_umask = FALSE)
  list(enriched = enriched, provenance = provenance, sha256 = derivative_sha256)
}

audit_aec_division_structure <- function(data) {
  required <- c("State", "DivisionID_2019", "DivisionNm_2019",
                "DivisionID_2022", "DivisionNm_2022")
  if (!all(required %in% names(data))) stop("Enriched division columns are incomplete", call. = FALSE)
  counts <- do.call(rbind, lapply(c("2019", "2022"), function(year) {
    id <- data[[paste0("DivisionID_", year)]]
    state_counts <- tapply(id, data$State, function(x) length(unique(x)))
    rbind(
      data.frame(clustering_year = year, State = names(state_counts),
                 unique_divisions = as.integer(state_counts), stringsAsFactors = FALSE),
      data.frame(clustering_year = year, State = "NATIONAL",
                 unique_divisions = length(unique(id)), stringsAsFactors = FALSE)
    )
  }))
  transitions <- aggregate(
    list(polling_place_rows = data$match_id),
    list(State = data$State, DivisionID_2019 = data$DivisionID_2019,
         DivisionNm_2019 = data$DivisionNm_2019,
         DivisionID_2022 = data$DivisionID_2022,
         DivisionNm_2022 = data$DivisionNm_2022),
    FUN = length
  )
  transitions <- transitions[order(transitions$State, transitions$DivisionID_2019,
                                   transitions$DivisionID_2022), , drop = FALSE]
  identical_rows <- data$DivisionID_2019 == data$DivisionID_2022
  summary <- data.frame(
    analytic_rows = nrow(data), rows_identical_division_id = sum(identical_rows),
    rows_changed_division_id = sum(!identical_rows),
    all_division_ids_identical = all(identical_rows),
    transition_cells = nrow(transitions), stringsAsFactors = FALSE
  )
  list(counts = counts, transitions = transitions, summary = summary)
}

aec_division_cluster_indices <- function(data, year) {
  division_column <- paste0("DivisionID_", year)
  if (!division_column %in% names(data)) stop("Unknown clustering year", call. = FALSE)
  state_rows <- split(seq_len(nrow(data)), data$State)
  lapply(state_rows, function(ii) split(ii, data[[division_column]][ii]))
}

aec_division_cluster_replicate <- function(task, data, cluster_indices) {
  with_rng_stream(task$seed, {
    warnings <- character()
    result <- tryCatch(withCallingHandlers({
      draws <- lapply(cluster_indices, function(divisions) {
        sample.int(length(divisions), length(divisions), replace = TRUE)
      })
      index <- unlist(Map(function(divisions, selected) {
        unlist(divisions[selected], use.names = FALSE)
      }, cluster_indices, draws), use.names = FALSE)
      expected_rows <- sum(unlist(Map(function(divisions, selected) {
        lengths(divisions)[selected]
      }, cluster_indices, draws), use.names = FALSE))
      X <- cbind(`2019` = data$X_2019[index], `2022` = data$X_2022[index])
      M <- cbind(`2019` = data$M_2019[index], `2022` = data$M_2022[index])
      fit <- fit_latent_binomial_cov(X, M, keep_influence = FALSE)
      rho <- fit$R_latent[1L, 2L]
      finite <- is.finite(rho) && abs(rho) < 1
      data.frame(
        clustering_year = as.integer(task$clustering_year),
        replicate = as.integer(task$replicate), stream_id = as.integer(task$stream_id),
        bootstrap_rows = length(index), expected_rows_from_complete_clusters = expected_rows,
        states = length(cluster_indices), divisions_drawn = sum(lengths(cluster_indices)),
        rho = rho, fisher = if (finite) atanh(rho) else NA_real_,
        projection_active = fit$projection_active, floor_active = fit$floor_active,
        finite = finite, complete_cluster_invariant = identical(length(index), expected_rows),
        success = finite && identical(length(index), expected_rows), failure_reason = "",
        warning_count = 0L, warning_messages = "", stringsAsFactors = FALSE
      )
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) data.frame(
      clustering_year = as.integer(task$clustering_year),
      replicate = as.integer(task$replicate), stream_id = as.integer(task$stream_id),
      bootstrap_rows = NA_integer_, expected_rows_from_complete_clusters = NA_integer_,
      states = length(cluster_indices), divisions_drawn = sum(lengths(cluster_indices)),
      rho = NA_real_, fisher = NA_real_, projection_active = NA, floor_active = NA,
      finite = FALSE, complete_cluster_invariant = FALSE, success = FALSE,
      failure_reason = conditionMessage(e), warning_count = 0L, warning_messages = "",
      stringsAsFactors = FALSE
    ))
    result$warning_count <- length(warnings)
    result$warning_messages <- paste(unique(warnings), collapse = " | ")
    result
  })
}

canonical_aec_division_results <- function(x) {
  if (!nrow(x)) return(x)
  x <- x[order(x$clustering_year, x$replicate), , drop = FALSE]
  rownames(x) <- NULL
  x
}

run_checkpointed_aec_division_tasks <- function(
    tasks, data, cluster_indices, checkpoint_path, final_path,
    cluster = NULL, checkpoint_every = 250L) {
  existing <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else data.frame()
  completed <- if (nrow(existing)) existing$replicate else integer()
  pending <- tasks[!vapply(tasks, function(task) task$replicate, integer(1)) %in% completed]
  if (!is.null(cluster)) {
    phase5c_worker_data <- data
    phase5c_worker_clusters <- cluster_indices
    parallel::clusterExport(cluster, c("phase5c_worker_data", "phase5c_worker_clusters"),
                            envir = environment())
  }
  starts <- if (length(pending)) seq.int(1L, length(pending), by = checkpoint_every) else integer()
  for (start in starts) {
    chunk <- pending[start:min(start + checkpoint_every - 1L, length(pending))]
    new <- if (is.null(cluster)) {
      lapply(chunk, aec_division_cluster_replicate, data = data,
             cluster_indices = cluster_indices)
    } else {
      parallel::parLapplyLB(cluster, chunk, function(task) {
        aec_division_cluster_replicate(task, phase5c_worker_data, phase5c_worker_clusters)
      })
    }
    existing <- if (nrow(existing)) rbind(existing, do.call(rbind, new)) else do.call(rbind, new)
    existing <- canonical_aec_division_results(existing)
    atomic_save_rds(existing, checkpoint_path)
  }
  complete <- nrow(existing) == length(tasks) && !anyDuplicated(existing$replicate) &&
    identical(existing$replicate, seq_len(length(tasks)))
  if (complete) atomic_save_rds(existing, final_path)
  list(results = existing, complete = complete, completed = nrow(existing), planned = length(tasks))
}

summarise_aec_division_bootstrap <- function(results, data, year, workers,
                                             master_seed = 20260909L) {
  valid <- results$success & results$finite & is.finite(results$fisher)
  z_limits <- quantile(results$fisher[valid], c(0.025, 0.975), type = 8, names = FALSE)
  interval <- tanh(z_limits)
  booth_interval <- c(0.908676730, 0.917612822)
  clusters <- aec_division_cluster_indices(data, year)
  state_counts <- lengths(clusters)
  data.frame(
    clustering_year = as.integer(year), states = length(clusters),
    divisions_by_state = paste(names(state_counts), state_counts, sep = "=", collapse = "|"),
    total_unique_divisions = length(unique(data[[paste0("DivisionID_", year)]])),
    planned_replicates = 4999L, valid_replicates = sum(valid),
    invalid_replicates = sum(!valid), nonfinite_count = sum(!results$finite),
    bootstrap_mean_correlation = mean(results$rho[valid]),
    bootstrap_median_correlation = median(results$rho[valid]),
    bootstrap_sd_fisher = sd(results$fisher[valid]),
    minimum_bootstrap_rows = min(results$bootstrap_rows[valid]),
    median_bootstrap_rows = median(results$bootstrap_rows[valid]),
    maximum_bootstrap_rows = max(results$bootstrap_rows[valid]),
    projection_activation_count = sum(results$projection_active %in% TRUE),
    floor_activation_count = sum(results$floor_active %in% TRUE),
    quantile_type = 8L, ci_lower = interval[1L], ci_upper = interval[2L],
    ci_width = diff(interval), booth_bootstrap_ci_width = diff(booth_interval),
    cluster_to_booth_ci_width_ratio = diff(interval) / diff(booth_interval),
    primary_point_inside_interval = interval[1L] <= 0.913227726 && interval[2L] >= 0.913227726,
    primary_point_reference = 0.913227726,
    analytic_ci_reference = "0.908574287|0.917654531",
    booth_bootstrap_ci_reference = "0.908676730|0.917612822",
    master_seed = as.integer(master_seed), worker_count = as.integer(workers),
    stringsAsFactors = FALSE
  )
}
