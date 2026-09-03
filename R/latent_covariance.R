fit_latent_binomial_cov <- function(X, M, eigen_floor = NULL,
                                    keep_influence = TRUE) {
  info <- validate_binomial_data(X, M)
  X <- .as_numeric_matrix(X, "X")
  M <- .as_numeric_matrix(M, "M")
  storage.mode(X) <- "double"
  storage.mode(M) <- "double"
  colnames(X) <- colnames(M) <- info$coordinate_names
  n <- info$n
  d <- info$d

  default_floor <- is.null(eigen_floor)
  if (default_floor) eigen_floor <- n^(-2)
  if (length(eigen_floor) != 1L || !is.finite(eigen_floor) || eigen_floor <= 0) {
    stop("eigen_floor must be one finite strictly positive number", call. = FALSE)
  }

  Y <- X / M
  Q <- Y * (1 - Y) / (M - 1)
  colnames(Y) <- colnames(Q) <- info$coordinate_names
  ybar <- colMeans(Y)
  qbar <- colMeans(Q)
  Yc <- sweep(Y, 2L, ybar, "-")
  S_Y <- crossprod(Yc) / (n - 1)
  Sigma_emp_n <- crossprod(Yc) / n - diag(qbar, nrow = d, ncol = d)
  Sigma_raw <- S_Y - diag(qbar, nrow = d, ncol = d)
  Sigma_raw <- symmetrise_matrix(Sigma_raw)

  matrix_names <- list(info$coordinate_names, info$coordinate_names)
  for (A_name in c("S_Y", "Sigma_emp_n", "Sigma_raw")) {
    A <- get(A_name)
    dimnames(A) <- matrix_names
    assign(A_name, A)
  }
  names(ybar) <- names(qbar) <- info$coordinate_names

  projection <- project_psd(Sigma_raw, eigen_floor = 0)
  lam <- projection$raw_eigenvalues
  V <- projection$eigenvectors
  stable_values <- pmax(lam, eigen_floor)
  Sigma_stable <- V %*% diag(stable_values, d, d) %*% t(V)
  Sigma_stable <- symmetrise_matrix(Sigma_stable)
  Sigma_psd <- projection$matrix
  dimnames(Sigma_psd) <- dimnames(Sigma_stable) <- matrix_names

  sd_latent <- sqrt(diag(Sigma_stable))
  R_latent <- Sigma_stable / outer(sd_latent, sd_latent)
  R_latent <- snap_unit_roundoff(R_latent)
  diag(R_latent) <- 1
  dimnames(R_latent) <- matrix_names

  Theta <- solve(Sigma_stable)
  Theta <- symmetrise_matrix(Theta)
  dimnames(Theta) <- matrix_names
  theta_sd <- sqrt(diag(Theta))
  Partial <- -Theta / outer(theta_sd, theta_sd)
  Partial <- snap_unit_roundoff(Partial)
  diag(Partial) <- 1
  dimnames(Partial) <- matrix_names

  reliability <- rep(NA_real_, d)
  positive_sy <- diag(S_Y) > 0
  reliability[positive_sy] <- pmin(1, pmax(0, diag(Sigma_psd)[positive_sy] /
                                            diag(S_Y)[positive_sy]))
  names(reliability) <- info$coordinate_names

  Psi_hat <- NULL
  influence_center_error <- NA_real_
  if (isTRUE(keep_influence)) {
    Psi_hat <- array(NA_real_, dim = c(n, d, d),
                     dimnames = list(NULL, info$coordinate_names, info$coordinate_names))
    for (i in seq_len(n)) {
      Psi_hat[i, , ] <- tcrossprod(Yc[i, ]) - diag(Q[i, ], d, d) - Sigma_emp_n
    }
    influence_center_error <- max(abs(apply(Psi_hat, c(2L, 3L), mean)))
  } else if (!identical(keep_influence, FALSE)) {
    stop("keep_influence must be TRUE or FALSE", call. = FALSE)
  }

  fit <- list(
    n = n, d = d, coordinate_names = info$coordinate_names,
    X = X, M = M, Y = Y, Q = Q, ybar = ybar, qbar = qbar,
    S_Y = S_Y, Sigma_emp_n = Sigma_emp_n, Sigma_raw = Sigma_raw,
    raw_eigenvalues = lam, min_eigen_raw = min(lam),
    Sigma_psd = Sigma_psd, Sigma_stable = Sigma_stable, lambda_n = eigen_floor,
    projection_active = any(lam < 0), floor_active = any(lam < eigen_floor),
    R_latent = R_latent, Theta = Theta, Partial = Partial,
    reliability = reliability, Psi_hat = Psi_hat,
    validation_summary = list(
      n = n, d = d, coordinate_names = info$coordinate_names,
      integer_tolerance = INTEGER_TOL,
      default_floor = default_floor,
      influence_center_error = influence_center_error,
      unweighted_rows = TRUE
    )
  )
  if (!isTRUE(keep_influence)) fit$Psi_hat <- NULL
  class(fit) <- "latent_binomial_cov"
  fit
}

print.latent_binomial_cov <- function(x, ...) {
  cat(sprintf("Latent binomial covariance fit: n=%d, d=%d\n", x$n, x$d))
  cat(sprintf("PSD projection active: %s; minimum raw eigenvalue: %.8g\n",
              x$projection_active, x$min_eigen_raw))
  cat("Stable covariance:\n")
  print(x$Sigma_stable)
  cat("Latent correlation:\n")
  print(x$R_latent)
  invisible(x)
}

summary.latent_binomial_cov <- function(object, ...) {
  out <- object[c(
    "n", "d", "coordinate_names", "ybar", "qbar", "S_Y", "Sigma_emp_n",
    "Sigma_raw", "raw_eigenvalues", "min_eigen_raw", "Sigma_psd",
    "Sigma_stable", "lambda_n", "projection_active", "floor_active",
    "R_latent", "Theta", "Partial", "reliability", "validation_summary"
  )]
  class(out) <- "summary.latent_binomial_cov"
  out
}

print.summary.latent_binomial_cov <- function(x, ...) {
  cat(sprintf("Latent binomial covariance summary: n=%d, d=%d\n", x$n, x$d))
  cat(sprintf("Projection active: %s; floor active: %s; lambda_n=%.8g\n",
              x$projection_active, x$floor_active, x$lambda_n))
  cat(sprintf("Minimum raw eigenvalue: %.8g\n", x$min_eigen_raw))
  cat("Reliability:\n")
  print(x$reliability)
  cat("Latent correlation:\n")
  print(x$R_latent)
  invisible(x)
}

coef.latent_binomial_cov <- function(object,
                                     type = c("covariance", "correlation", "partial"), ...) {
  type <- match.arg(type)
  switch(type,
         covariance = object$Sigma_stable,
         correlation = object$R_latent,
         partial = object$Partial)
}
