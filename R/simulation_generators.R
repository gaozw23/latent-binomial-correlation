read_pairwise_scenarios <- function(path = file.path(
                                      getOption("latent_binomial_project_root", "."),
                                      "config", "pairwise_scenarios.csv")) {
  scenarios <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("scenario_id", "latent_sign", "a", "b", "delta",
                "denominator_id", "description")
  if (!identical(names(scenarios), required)) {
    stop("Pairwise scenario configuration has an invalid schema", call. = FALSE)
  }
  scenarios
}

.scenario_row <- function(scenario) {
  if (is.character(scenario) && length(scenario) == 1L) {
    scenarios <- read_pairwise_scenarios()
    hit <- scenarios[scenarios$scenario_id == scenario, , drop = FALSE]
    if (nrow(hit) != 1L) stop("scenario_id must identify exactly one configured scenario", call. = FALSE)
    return(hit)
  }
  if (is.data.frame(scenario) && nrow(scenario) == 1L) return(scenario)
  stop("scenario must be one configured scenario_id or one scenario data-frame row", call. = FALSE)
}

pairwise_truth <- function(scenario) {
  s <- .scenario_row(scenario)
  mu1 <- s$a / (s$a + s$b)
  if (s$latent_sign == "positive") {
    mu2 <- mu1
    sign <- 1
  } else if (s$latent_sign == "negative") {
    mu2 <- (1 - s$delta) * (1 - mu1) + s$delta * mu1
    sign <- -1
  } else {
    stop("latent_sign must be 'positive' or 'negative'", call. = FALSE)
  }
  rho <- sign * (1 - s$delta) / sqrt((1 - s$delta)^2 + s$delta^2)
  variance <- s$a * s$b / ((s$a + s$b)^2 * (s$a + s$b + 1))
  Sigma <- variance * matrix(c(
    1,
    sign * (1 - s$delta),
    sign * (1 - s$delta),
    (1 - s$delta)^2 + s$delta^2
  ), 2L, 2L)
  list(mu = c(mu1, mu2), Sigma = Sigma, rho = rho)
}

generate_pairwise_latent <- function(n, scenario) {
  s <- .scenario_row(scenario)
  if (length(n) != 1L || n < 1L || n != as.integer(n)) stop("n must be a positive integer", call. = FALSE)
  U <- rbeta(n, s$a, s$b)
  V <- rbeta(n, s$a, s$b)
  if (s$latent_sign == "positive") {
    P1 <- U
    P2 <- (1 - s$delta) * U + s$delta * V
  } else if (s$latent_sign == "negative") {
    P1 <- U
    P2 <- (1 - s$delta) * (1 - U) + s$delta * V
  } else {
    stop("latent_sign must be 'positive' or 'negative'", call. = FALSE)
  }
  cbind(P1 = P1, P2 = P2)
}

generate_pairwise_data <- function(n, scenario) {
  s <- .scenario_row(scenario)
  P <- generate_pairwise_latent(n, s)
  truth <- pairwise_truth(s)
  P1 <- P[, 1L]
  P2 <- P[, 2L]
  if (s$denominator_id == "D1") {
    N1 <- rpois(n, 18)
    N2 <- rpois(n, 18)
    M1 <- 2L + N1
    M2 <- 2L + N2
    denominator_components <- list(N1 = N1, N2 = N2)
  } else if (s$denominator_id == "D2") {
    N1 <- rpois(n, 2)
    N2 <- rpois(n, 20)
    M1 <- 2L + N1
    M2 <- 2L + N2
    denominator_components <- list(N1 = N1, N2 = N2)
  } else if (s$denominator_id == "D3") {
    H <- ((P1 - truth$mu[1L]) + (P2 - truth$mu[2L])) / 2
    C <- rpois(n, 3 * exp(1.5 * H))
    E1 <- rpois(n, 4 * exp(2 * (P1 - truth$mu[1L])))
    E2 <- rpois(n, 12 * exp(-1.5 * (P2 - truth$mu[2L])))
    M1 <- 2L + C + E1
    M2 <- 2L + C + E2
    denominator_components <- list(C = C, E1 = E1, E2 = E2,
                                   shared_C_used_once = TRUE)
  } else {
    stop("Unknown denominator_id", call. = FALSE)
  }
  X1 <- rbinom(n, M1, P1)
  X2 <- rbinom(n, M2, P2)
  list(
    X = cbind(V1 = X1, V2 = X2),
    M = cbind(V1 = M1, V2 = M2),
    P = P,
    truth = truth,
    scenario = s,
    denominator_components = denominator_components
  )
}

matrix_weight_matrix <- function() {
  rbind(
    c(0.25, 0.35, 0.00, 0.40, 0.00, 0.00, 0.00, 0.00),
    c(0.25, 0.35, 0.00, 0.00, 0.40, 0.00, 0.00, 0.00),
    c(0.25, 0.20, 0.15, 0.00, 0.00, 0.40, 0.00, 0.00),
    c(0.25, 0.00, 0.35, 0.00, 0.00, 0.00, 0.40, 0.00),
    c(0.25, 0.00, 0.35, 0.00, 0.00, 0.00, 0.00, 0.40)
  )
}

matrix_truth <- function() {
  W <- matrix_weight_matrix()
  if (any(W < 0) || any(rowSums(W) != 1)) stop("W must be nonnegative with exact unit row sums", call. = FALSE)
  mu_B <- 2 / 7
  sigma2_B <- 2 * 5 / ((2 + 5)^2 * (2 + 5 + 1))
  Sigma <- sigma2_B * W %*% t(W)
  Theta <- solve(Sigma)
  Partial <- -Theta / sqrt(outer(diag(Theta), diag(Theta)))
  diag(Partial) <- 1
  list(W = W, mu = rep(mu_B, 5L), sigma2_B = sigma2_B,
       Sigma = Sigma, R = cov2cor(Sigma), Theta = Theta, Partial = Partial)
}

generate_matrix_latent <- function(n) {
  if (length(n) != 1L || n < 1L || n != as.integer(n)) stop("n must be a positive integer", call. = FALSE)
  F <- matrix(NA_real_, nrow = n, ncol = 8L)
  for (j in seq_len(8L)) F[, j] <- rbeta(n, 2, 5)
  colnames(F) <- c("G", "A", "B", "E1", "E2", "E3", "E4", "E5")
  P <- F %*% t(matrix_weight_matrix())
  colnames(P) <- paste0("V", seq_len(5L))
  list(F = F, P = P)
}

generate_matrix_data <- function(n) {
  latent <- generate_matrix_latent(n)
  P <- latent$P
  mu_B <- 2 / 7
  lambda <- c(3, 5, 8, 12, 20)
  gamma <- c(1.5, -1.25, 1, 0.75, -0.5)
  C <- rpois(n, 4 * exp(1.2 * (rowMeans(P) - mu_B)))
  EM <- M <- X <- matrix(NA_real_, nrow = n, ncol = 5L)
  for (j in seq_len(5L)) {
    EM[, j] <- rpois(n, lambda[j] * exp(gamma[j] * (P[, j] - mu_B)))
    M[, j] <- 2L + C + EM[, j]
  }
  for (j in seq_len(5L)) X[, j] <- rbinom(n, M[, j], P[, j])
  colnames(X) <- colnames(M) <- colnames(P)
  list(X = X, M = M, P = P, F = latent$F, EM = EM, C = C, truth = matrix_truth())
}
