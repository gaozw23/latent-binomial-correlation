PROJECT_ROOT <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
options(latent_binomial_project_root = PROJECT_ROOT)
Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME"))
source(file.path(PROJECT_ROOT, "R", "constants.R"), local = FALSE)
source(file.path(PROJECT_ROOT, "R", "utilities.R"), local = FALSE)
files <- sort(list.files(file.path(PROJECT_ROOT, "R"), pattern = "[.]R$", full.names = TRUE))
invisible(lapply(files, sys.source, envir = .GlobalEnv))

core_test_streams <- make_streams(11L, MASTER_SEEDS[["core_unit_tests"]])

record_phase1_diagnostic <- function(test_id, value) {
  audit_dir <- file.path(PROJECT_ROOT, "results", "audits")
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, file.path(audit_dir, paste0(test_id, "_diagnostic.rds")))
  invisible(value)
}

weighted_covariance_functional <- function(Y, Q, weights) {
  weights <- weights / sum(weights)
  mu <- colSums(Y * weights)
  Yc <- sweep(Y, 2L, mu, "-")
  Sigma <- crossprod(Yc, Yc * weights) - diag(colSums(Q * weights), ncol(Y))
  list(mu = mu, Sigma = symmetrise_matrix(Sigma))
}
