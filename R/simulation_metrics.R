raw_covariance_estimate <- function(X, M) {
  validate_binomial_data(X, M)
  X <- .as_numeric_matrix(X, "X")
  M <- .as_numeric_matrix(M, "M")
  Y <- X / M
  Q <- Y * (1 - Y) / (M - 1)
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  symmetrise_matrix(crossprod(Yc) / (nrow(Y) - 1) - diag(colMeans(Q), ncol(Y)))
}

frobenius_projection_check <- function(Sigma_raw, Sigma_psd, Sigma_target,
                                       tolerance = PROJECTION_LOSS_TOL) {
  raw_loss <- norm(Sigma_raw - Sigma_target, "F")
  psd_loss <- norm(Sigma_psd - Sigma_target, "F")
  allowed <- raw_loss + tolerance * (1 + raw_loss)
  list(pass = psd_loss <= allowed, raw_loss = raw_loss,
       psd_loss = psd_loss, margin = allowed - psd_loss)
}
