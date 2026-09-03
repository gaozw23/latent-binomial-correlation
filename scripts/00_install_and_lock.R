options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

if (!file.exists("renv/activate.R")) {
  renv::init(bare = TRUE, restart = FALSE)
}

required <- c(
  "Matrix", "data.table", "ggplot2", "patchwork", "testthat", "jsonlite",
  "digest", "curl", "xml2", "rvest", "stringi", "future", "future.apply",
  "withr"
)
renv::hydrate(packages = required, prompt = FALSE)
renv::snapshot(packages = c("renv", required), prompt = FALSE)
renv::restore(prompt = FALSE)

cat("renv installation, snapshot, and restoration completed successfully\n")
