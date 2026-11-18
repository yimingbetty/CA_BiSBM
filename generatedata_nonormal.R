## Additional data generators for the four-case experiment.
## This file does not modify generatedata.R; it sources it and adds new objects.

source("generatedata.R")

## Case 2: unbalanced Bernoulli block probabilities.
B_unbalanced <- matrix(c(
  0.160, 0.085, 0.055, 0.035, 0.020,
  0.075, 0.125, 0.095, 0.060, 0.040,
  0.020, 0.040, 0.075, 0.145, 0.180
), nrow = 3, byrow = TRUE)

## Heavy-tailed non-Gaussian covariates.
generate_nonnormal_covariates <- function(
  n, K, labels, p, M, sigma = 1, df = 3
) {
  if (df <= 2) stop("df must be greater than 2.")
  if (!all(dim(M) == c(K, p))) stop("M has incompatible dimensions.")

  X <- matrix(NA_real_, n, p)
  t_scale <- sqrt(df / (df - 2))

  for (i in seq_len(n)) {
    noise <- rt(p, df = df) / t_scale
    X[i, ] <- M[labels[i], ] + sigma * noise
  }

  list(
    X = X,
    M = M,
    Sigma_list = lapply(seq_len(K), function(k) diag(sigma^2, p))
  )
}

generate_nonnormal_ca_bisbm <- function(
  n_r, n_c, K_r, K_c, B, p_r, p_c,
  M_r = NULL, M_c = NULL,
  sigma_r = 1, sigma_c = 1,
  df_r = 3, df_c = 3,
  fixed = NULL
) {
  dat <- generate_ca_bisbm(
    n_r = n_r, n_c = n_c,
    K_r = K_r, K_c = K_c,
    B = B, p_r = p_r, p_c = p_c,
    M_r = M_r, M_c = M_c,
    sigma_r = sigma_r, sigma_c = sigma_c,
    fixed = fixed
  )

  M_r_use <- if (!is.null(fixed) && !is.null(fixed$M_r)) fixed$M_r else M_r
  M_c_use <- if (!is.null(fixed) && !is.null(fixed$M_c)) fixed$M_c else M_c

  if (is.null(M_r_use) || is.null(M_c_use)) {
    stop("Case 3 requires M_r and M_c.")
  }

  X_res <- generate_nonnormal_covariates(
    n_r, K_r, dat$c_r, p_r, M_r_use, sigma_r, df_r
  )
  Y_res <- generate_nonnormal_covariates(
    n_c, K_c, dat$c_c, p_c, M_c_use, sigma_c, df_c
  )

  dat$X <- X_res$X
  dat$Y <- Y_res$X
  dat$M_r <- X_res$M
  dat$M_c <- Y_res$M
  dat$Sigma_r_list <- X_res$Sigma_list
  dat$Sigma_c_list <- Y_res$Sigma_list
  dat$covariate_distribution <- "student_t"
  dat$df_r <- df_r
  dat$df_c <- df_c
  dat
}

generate_connected_nonnormal_ca_bisbm <- function(..., max_tries = 1000) {
  args <- list(...)
  for (attempt in seq_len(max_tries)) {
    dat <- do.call(generate_nonnormal_ca_bisbm, args)
    if (all(rowSums(dat$A) > 0) && all(colSums(dat$A) > 0)) {
      dat$n_tries <- attempt
      return(dat)
    }
  }
  stop(sprintf(
    "Unable to generate a connected non-Gaussian network in %d attempts.",
    max_tries
  ))
}
