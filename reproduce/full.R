if (!identical(Sys.getenv("FULL_REPRODUCTION"), "YES")) {
  stop("Set FULL_REPRODUCTION=YES after reading the README and scripts", call. = FALSE)
}

steps <- c(
  "scripts/01_run_unit_tests.R",
  "scripts/16_generate_phase3_rng_streams.R",
  "scripts/04_pairwise_full.R",
  "scripts/05_bootstrap_interval_study.R",
  "scripts/18_generate_phase4_rng_streams.R",
  "scripts/07_matrix_full.R",
  "scripts/08_download_aec.R",
  "scripts/09_prepare_aec.R",
  "scripts/21_adjudicate_aec_outliers.R",
  "scripts/22_freeze_aec_matches.R",
  "scripts/10_analyse_aec.R",
  "scripts/11_controlled_thinning.R",
  "scripts/24_enrich_aec_divisions.R",
  "scripts/25_aec_division_cluster_bootstrap.R",
  "reproduce/quick.R"
)

missing <- steps[!file.exists(steps)]
if (length(missing)) stop("Missing reproduction files: ", paste(missing, collapse = ", "), call. = FALSE)
for (step in steps) {
  message("Running ", step)
  source(step, local = new.env(parent = globalenv()))
}
