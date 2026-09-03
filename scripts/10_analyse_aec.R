start_time <- Sys.time()
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

dir.create("results/summaries", recursive = TRUE, showWarnings = FALSE)
dir.create("results/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("results/audits", recursive = TRUE, showWarnings = FALSE)
dir.create("config/rng_streams", recursive = TRUE, showWarnings = FALSE)

data <- read_frozen_aec()
data$harmonic_M <- 2 / (1 / data$M_2019 + 1 / data$M_2022)
small_threshold <- as.numeric(quantile(data$harmonic_M, 0.25, type = 7))

primary <- analyse_aec(data, "primary", small_place_threshold = small_threshold)
primary_latent <- primary$summary$latent_correlation
exact <- analyse_aec(data[data$match_route == "EXACT_ID", , drop = FALSE],
                     "exact_id_only", primary_latent, small_threshold)
m20 <- analyse_aec(data[data$M_2019 >= 20 & data$M_2022 >= 20, , drop = FALSE],
                   "M_at_least_20", primary_latent, small_threshold)
small <- analyse_aec(data[data$harmonic_M <= small_threshold, , drop = FALSE],
                     "small_place_lower_quartile_ties_included", primary_latent, small_threshold)

exclusions <- read.csv("results/derived/aec_match_exclusions.csv", stringsAsFactors = FALSE)
primary$summary$M_lt_2_exclusions <- sum(exclusions$ExclusionReason == "M_LT_2_IN_AT_LEAST_ONE_YEAR")
write.csv(primary$summary, "results/summaries/aec_primary_analysis.csv", row.names = FALSE, quote = TRUE)
write.csv(exact$summary, "results/summaries/aec_exact_id_sensitivity.csv", row.names = FALSE, quote = TRUE)
write.csv(m20$summary, "results/summaries/aec_m20_sensitivity.csv", row.names = FALSE, quote = TRUE)
write.csv(small$summary, "results/summaries/aec_small_place_sensitivity.csv", row.names = FALSE, quote = TRUE)
saveRDS(list(primary = primary, exact_id = exact, m20 = m20, small_place = small,
             small_place_threshold = small_threshold),
        "results/raw/aec_analysis_fits.rds", compress = "xz")

bootstrap_streams <- make_streams(4999L, MASTER_SEEDS[["aec_bootstrap"]])
bootstrap_map <- data.frame(
  stream_id = seq_len(4999L), replicate = seq_len(4999L),
  study_component = "aec_state_stratified_bootstrap", stringsAsFactors = FALSE
)
saveRDS(bootstrap_streams, "config/rng_streams/aec_bootstrap_streams.rds", compress = "xz")
write.csv(bootstrap_map, "config/rng_streams/aec_bootstrap_stream_map.csv", row.names = FALSE, quote = TRUE)
bootstrap <- run_aec_state_bootstrap(data, bootstrap_streams)
bootstrap$summary$primary_point_estimate <- primary$summary$latent_correlation
bootstrap$summary$analytic_ci_lower <- primary$summary$analytic_ci_lower
bootstrap$summary$analytic_ci_upper <- primary$summary$analytic_ci_upper
saveRDS(bootstrap$replicates, "results/raw/aec_state_bootstrap.rds", compress = "xz")
write.csv(bootstrap$summary, "results/summaries/aec_state_bootstrap.csv", row.names = FALSE, quote = TRUE)

metadata <- list(
  started_utc = format(start_time, tz = "UTC", usetz = TRUE),
  ended_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
  frozen_rows = nrow(data), small_place_threshold = small_threshold,
  bootstrap_replicates = nrow(bootstrap$replicates)
)
saveRDS(metadata, "results/summaries/aec_analysis_execution_metadata.rds", compress = "xz")
cat(sprintf("AEC primary and sensitivities complete: n=%d, rho_naive=%.9f, rho_latent=%.9f\n",
            primary$summary$n, primary$summary$naive_correlation,
            primary$summary$latent_correlation))
cat(sprintf("State-stratified bootstrap: %d valid / %d; CI [%.9f, %.9f]\n",
            bootstrap$summary$valid_replicates, bootstrap$summary$bootstrap_replicates,
            bootstrap$summary$bootstrap_ci_lower, bootstrap$summary$bootstrap_ci_upper))
