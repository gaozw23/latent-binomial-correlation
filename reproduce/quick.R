if (!file.exists("scripts/27_build_phase6b_artifacts.R")) {
  stop("Run this script from the repository root", call. = FALSE)
}

old_output <- Sys.getenv("PHASE6B_OUTPUT_DIR", unset = NA_character_)
on.exit({
  if (is.na(old_output)) Sys.unsetenv("PHASE6B_OUTPUT_DIR") else
    Sys.setenv(PHASE6B_OUTPUT_DIR = old_output)
}, add = TRUE)

Sys.setenv(PHASE6B_OUTPUT_DIR = file.path("reproduce", "output"))
source("scripts/27_build_phase6b_artifacts.R")

cat("Quick reproduction complete: reproduce/output\n")

