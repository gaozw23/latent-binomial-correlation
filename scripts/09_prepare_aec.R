Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

dir.create("results/audits", recursive = TRUE, showWarnings = FALSE)
inventory <- parse_aec_inventory("results/manifests/aec_acquisition_manifest.csv")
audit <- audit_aec_inventory(inventory)
matching <- audit_aec_matching_feasibility(audit)

write.csv(audit$schema_inventory, "results/audits/aec_schema_inventory.csv",
          row.names = FALSE, quote = TRUE)
write.csv(audit$party_inventory, "results/audits/aec_party_inventory.csv",
          row.names = FALSE, quote = TRUE)
write.csv(audit$polling_place_type_inventory,
          "results/audits/aec_polling_place_type_inventory.csv",
          row.names = FALSE, quote = TRUE)
write.csv(audit$metadata_summary, "results/audits/aec_metadata_summary.csv",
          row.names = FALSE, quote = TRUE)
write.csv(audit$denominator_audit, "results/audits/aec_denominator_audit.csv",
          row.names = FALSE, quote = TRUE)
write.csv(audit$low_denominator_groups, "results/audits/aec_denominator_low_groups.csv",
          row.names = FALSE, quote = TRUE)
all_anomalies <- rbind(audit$anomaly_inventory, matching$anomaly_inventory)
write.csv(all_anomalies, "results/audits/aec_anomaly_inventory.csv",
          row.names = FALSE, quote = TRUE)
write.csv(matching$feasibility, "results/audits/aec_matching_feasibility.csv",
          row.names = FALSE, quote = TRUE)
write.csv(matching$exact_candidates, "results/audits/aec_exact_id_candidates.csv",
          row.names = FALSE, quote = TRUE)
write.csv(matching$secondary_candidates, "results/audits/aec_secondary_match_candidates.csv",
          row.names = FALSE, quote = TRUE)
write.csv(matching$ambiguous_names, "results/audits/aec_ambiguous_name_candidates.csv",
          row.names = FALSE, quote = TRUE)
write.csv(matching$distance_distribution,
          "results/audits/aec_secondary_distance_distribution.csv",
          row.names = FALSE, quote = TRUE)
write.csv(matching$near_threshold, "results/audits/aec_secondary_near_2km.csv",
          row.names = FALSE, quote = TRUE)

cat("AEC schema and matching-feasibility audit complete\n")
print(audit$metadata_summary)
print(matching$feasibility[matching$feasibility$state_ab == "NATIONAL", ])
