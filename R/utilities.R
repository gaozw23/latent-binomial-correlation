make_streams <- function(N, master_seed) {
  if (length(N) != 1L || is.na(N) || N < 1 || N != as.integer(N)) {
    stop("N must be one positive integer", call. = FALSE)
  }
  if (length(master_seed) != 1L || is.na(master_seed)) {
    stop("master_seed must be one non-missing scalar", call. = FALSE)
  }
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(master_seed))
  out <- vector("list", as.integer(N))
  out[[1L]] <- .Random.seed
  if (N >= 2L) {
    for (i in 2L:N) out[[i]] <- parallel::nextRNGStream(out[[i - 1L]])
  }
  out
}

with_rng_stream <- function(stream, code) {
  if (!is.numeric(stream) || length(stream) < 2L) {
    stop("stream must be a complete .Random.seed vector", call. = FALSE)
  }
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  assign(".Random.seed", stream, envir = .GlobalEnv)
  eval.parent(substitute(code))
}

source_project_r <- function(envir = parent.frame()) {
  files <- sort(list.files("R", pattern = "[.]R$", full.names = TRUE))
  invisible(lapply(files, sys.source, envir = envir))
}

sha256_files <- function(paths) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for SHA-256 manifests", call. = FALSE)
  }
  paths <- sort(normalizePath(paths, winslash = "/", mustWork = TRUE))
  data.frame(
    path = paths,
    sha256 = vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE),
    stringsAsFactors = FALSE
  )
}

snap_unit_roundoff <- function(x, tolerance = ROUND_OFF_SNAP_TOL) {
  for (target in c(-1, 0, 1)) {
    hit <- abs(x - target) <= tolerance
    x[hit] <- target
  }
  x
}

symmetrise_matrix <- function(A) (A + t(A)) / 2
