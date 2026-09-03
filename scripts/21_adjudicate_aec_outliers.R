Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

inventory <- parse_aec_inventory("results/manifests/aec_acquisition_manifest.csv")
audit <- audit_aec_inventory(inventory)
disposition_path <- "config/aec_exact_id_outlier_dispositions.csv"
adjudication <- build_aec_exact_id_outlier_adjudication(
  audit, if (file.exists(disposition_path)) disposition_path else NULL
)
dir.create("results/audits", recursive = TRUE, showWarnings = FALSE)
write.csv(adjudication, "results/audits/aec_exact_id_outlier_adjudication.csv",
          row.names = FALSE, quote = TRUE)
cat("Exact-ID adjudication rows:", nrow(adjudication), "\n")
print(table(adjudication$ProposedDisposition))
