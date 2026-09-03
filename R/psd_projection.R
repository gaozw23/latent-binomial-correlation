project_psd <- function(A, eigen_floor = 0, symmetry_tol = SYMMETRY_TOL) {
  if (!is.matrix(A) || !is.numeric(A) || is.complex(A) || nrow(A) != ncol(A)) {
    stop("A must be a real numeric square matrix", call. = FALSE)
  }
  if (nrow(A) < 1L || any(!is.finite(A))) {
    stop("A must be non-empty and contain only finite values", call. = FALSE)
  }
  if (length(eigen_floor) != 1L || !is.finite(eigen_floor) || eigen_floor < 0) {
    stop("eigen_floor must be one finite nonnegative number", call. = FALSE)
  }
  if (length(symmetry_tol) != 1L || !is.finite(symmetry_tol) || symmetry_tol < 0) {
    stop("symmetry_tol must be one finite nonnegative number", call. = FALSE)
  }

  asymmetry <- norm(A - t(A), "F")
  allowed <- symmetry_tol * (1 + norm(A, "F"))
  if (asymmetry > allowed) {
    stop(sprintf("A asymmetry %.17g exceeds allowed tolerance %.17g", asymmetry, allowed),
         call. = FALSE)
  }

  A_sym <- symmetrise_matrix(A)
  e <- eigen(A_sym, symmetric = TRUE)
  adjusted <- pmax(e$values, eigen_floor)
  d <- nrow(A)
  projected <- e$vectors %*% diag(adjusted, nrow = d, ncol = d) %*% t(e$vectors)
  projected <- symmetrise_matrix(projected)
  dimnames(projected) <- dimnames(A)

  list(
    matrix = projected,
    raw_eigenvalues = e$values,
    adjusted_eigenvalues = adjusted,
    eigenvectors = e$vectors,
    min_eigen_raw = min(e$values),
    projection_active = isTRUE(eigen_floor == 0) && any(e$values < 0),
    floor_active = any(e$values < eigen_floor),
    frobenius_adjustment = norm(projected - A_sym, "F"),
    original_asymmetry = asymmetry,
    symmetrised_input = A_sym
  )
}
