.first_matrix_failure <- function(mask) {
  hits <- which(mask, arr.ind = TRUE)
  hits <- hits[order(hits[, 1L], hits[, 2L]), , drop = FALSE]
  c(row = hits[1L, 1L], column = hits[1L, 2L])
}

.validation_error <- function(object, value, row, column, rule, coordinate_names = NULL) {
  column_label <- if (is.null(coordinate_names)) {
    as.character(column)
  } else {
    sprintf("%d ('%s')", column, coordinate_names[column])
  }
  stop(
    sprintf("%s validation failed at row %d, column %s: value=%s; rule=%s",
            object, row, column_label, format(value, digits = 17), rule),
    call. = FALSE
  )
}

.as_numeric_matrix <- function(x, object) {
  if (is.data.frame(x)) {
    ok <- vapply(x, function(col) is.numeric(col) || is.integer(col), logical(1))
    if (!all(ok)) {
      j <- which(!ok)[1L]
      stop(sprintf("%s column %d ('%s') must be numeric or integer",
                   object, j, names(x)[j]), call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !(is.numeric(x) || is.integer(x)) || is.complex(x)) {
    stop(sprintf("%s must be a numeric matrix or an all-numeric data frame", object),
         call. = FALSE)
  }
  x
}

validate_binomial_data <- function(X, M, require_integer = TRUE) {
  if (!is.logical(require_integer) || length(require_integer) != 1L || is.na(require_integer)) {
    stop("require_integer must be TRUE or FALSE", call. = FALSE)
  }
  X <- .as_numeric_matrix(X, "X")
  M <- .as_numeric_matrix(M, "M")

  if (!identical(dim(X), dim(M))) {
    stop(sprintf("X and M must have identical dimensions; got %s and %s",
                 paste(dim(X), collapse = "x"), paste(dim(M), collapse = "x")),
         call. = FALSE)
  }
  n <- nrow(X)
  d <- ncol(X)
  if (n < 2L) stop(sprintf("X and M require nrow >= 2; got %d", n), call. = FALSE)
  if (d < 2L) stop(sprintf("X and M require ncol >= 2; got %d", d), call. = FALSE)

  x_names <- colnames(X)
  m_names <- colnames(M)
  if (!is.null(x_names) && !is.null(m_names) && !identical(x_names, m_names)) {
    stop("X and M column names must agree exactly when both are supplied", call. = FALSE)
  }
  coordinate_names <- if (!is.null(x_names)) x_names else if (!is.null(m_names)) m_names else paste0("V", seq_len(d))
  if (anyNA(coordinate_names) || any(!nzchar(coordinate_names)) || anyDuplicated(coordinate_names)) {
    stop("Coordinate names must be non-missing, non-empty, and unique", call. = FALSE)
  }

  checks <- list(
    list(object = "X", data = X, mask = !is.finite(X), rule = "all X entries must be finite"),
    list(object = "M", data = M, mask = !is.finite(M), rule = "all M entries must be finite")
  )
  if (isTRUE(require_integer)) {
    checks <- c(checks, list(
      list(object = "X", data = X, mask = is.finite(X) & abs(X - round(X)) > INTEGER_TOL,
           rule = sprintf("X must be integer-valued within %.0e", INTEGER_TOL)),
      list(object = "M", data = M, mask = is.finite(M) & abs(M - round(M)) > INTEGER_TOL,
           rule = sprintf("M must be integer-valued within %.0e", INTEGER_TOL))
    ))
  }
  checks <- c(checks, list(
    list(object = "M", data = M, mask = is.finite(M) & M < 2,
         rule = "M must be at least 2"),
    list(object = "X", data = X, mask = is.finite(X) & X < 0,
         rule = "X must be nonnegative"),
    list(object = "X", data = X, mask = is.finite(X) & is.finite(M) & X > M,
         rule = "X must not exceed M")
  ))

  for (check in checks) {
    if (any(check$mask)) {
      at <- .first_matrix_failure(check$mask)
      .validation_error(check$object, check$data[at[1L], at[2L]], at[1L], at[2L],
                        check$rule, coordinate_names)
    }
  }

  invisible(list(n = n, d = d, coordinate_names = coordinate_names))
}
