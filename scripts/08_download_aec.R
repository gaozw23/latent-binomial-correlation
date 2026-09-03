Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source("R/constants.R", local = .GlobalEnv)
source("R/utilities.R", local = .GlobalEnv)
source_project_r(.GlobalEnv)

dir.create("results/audits", recursive = TRUE, showWarnings = FALSE)
dir.create("results/manifests", recursive = TRUE, showWarnings = FALSE)

discovered <- discover_aec_links("config/aec_manifest.csv")
write.csv(discovered, "results/audits/aec_discovered_link_inventory.csv",
          row.names = FALSE, quote = TRUE)
manifest <- download_aec_files(discovered, "data/raw/aec")
write.csv(manifest, "results/manifests/aec_acquisition_manifest.csv",
          row.names = FALSE, quote = TRUE)

cat("AEC acquisition complete:", nrow(manifest), "/18 files\n")
print(table(manifest$election_year, manifest$file_role))
