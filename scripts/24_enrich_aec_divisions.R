Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

result <- create_aec_division_enrichment()
structure <- audit_aec_division_structure(result$enriched)
write.csv(structure$counts, "results/audits/aec_division_counts.csv",
          row.names = FALSE, quote = TRUE)
write.csv(structure$transitions, "results/audits/aec_division_transitions.csv",
          row.names = FALSE, quote = TRUE)
write.csv(structure$summary, "results/audits/aec_division_structure_summary.csv",
          row.names = FALSE, quote = TRUE)
cat(sprintf("Enriched rows: %d; SHA-256: %s\n", nrow(result$enriched), result$sha256))
cat(sprintf("Division IDs identical rows: %d; changed rows: %d\n",
            structure$summary$rows_identical_division_id,
            structure$summary$rows_changed_division_id))
