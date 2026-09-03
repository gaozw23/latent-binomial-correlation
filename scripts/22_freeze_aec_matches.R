Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

inventory <- parse_aec_inventory("results/manifests/aec_acquisition_manifest.csv")
audit <- audit_aec_inventory(inventory)
adjudication <- read.csv("results/audits/aec_exact_id_outlier_adjudication.csv",
                         stringsAsFactors = FALSE, colClasses = "character")
if (any(adjudication$ProposedDisposition == "UNRESOLVED")) {
  stop("Gate A failed: unresolved exact-ID adjudication", call. = FALSE)
}
result <- freeze_aec_matched_table(audit, adjudication)
cat("Frozen eligible matches:", nrow(result$frozen), "\n")
cat("Candidate matches:", result$candidate_matches,
    "exact:", result$exact_candidates, "secondary:", result$secondary_candidates, "\n")
cat("M<2 exclusions:", result$denominator_exclusions, "\n")
