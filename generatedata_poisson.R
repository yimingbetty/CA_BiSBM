

source("generatedata.R")

generate_poisson_nodc_bisbm <- function(
  n_r, n_c, K_r, K_c, B, rho = 0.8,
  p_r, p_c,
  M_r = NULL, M_c = NULL,
  Sigma_r_list = NULL, Sigma_c_list = NULL,
  sigma_r = 1, sigma_c = 1,
  row_cluster_groups = NULL, col_cluster_groups = NULL,
  fixed = NULL
) {
  if (!is.matrix(B) || !all(dim(B) == c(K_r, K_c))) {
    stop("B 必须是 K_r x K_c 的矩阵。")
  }
  if (any(B < 0) || !is.finite(rho) || rho <= 0) {
    stop("B 的元素必须非负，rho 必须为正数。")
  }

  if (!is.null(fixed) && !is.null(fixed$c_r) && !is.null(fixed$c_c)) {
    c_r <- fixed$c_r
    c_c <- fixed$c_c
  } else {
    c_r <- sample.int(K_r, n_r, replace = TRUE)
    c_c <- sample.int(K_c, n_c, replace = TRUE)
  }

  Z_r <- matrix(0, nrow = n_r, ncol = K_r)
  Z_c <- matrix(0, nrow = n_c, ncol = K_c)
  Z_r[cbind(seq_len(n_r), c_r)] <- 1
  Z_c[cbind(seq_len(n_c), c_c)] <- 1

 
  Lambda <- rho * (Z_r %*% B %*% t(Z_c))
  A <- matrix(rpois(n_r * n_c, lambda = as.vector(Lambda)), nrow = n_r, ncol = n_c)

  M_r_use <- if (!is.null(fixed) && !is.null(fixed$M_r)) fixed$M_r else M_r
  M_c_use <- if (!is.null(fixed) && !is.null(fixed$M_c)) fixed$M_c else M_c
  Sigma_r_use <- if (!is.null(fixed) && !is.null(fixed$Sigma_r_list)) {
    fixed$Sigma_r_list
  } else {
    Sigma_r_list
  }
  Sigma_c_use <- if (!is.null(fixed) && !is.null(fixed$Sigma_c_list)) {
    fixed$Sigma_c_list
  } else {
    Sigma_c_list
  }

  X_res <- generate_covariates_with_groups(
    n_r, K_r, c_r, p_r, M = M_r_use, Sigma_list = Sigma_r_use,
    sigma = sigma_r, cluster_groups = row_cluster_groups
  )
  Y_res <- generate_covariates_with_groups(
    n_c, K_c, c_c, p_c, M = M_c_use, Sigma_list = Sigma_c_use,
    sigma = sigma_c, cluster_groups = col_cluster_groups
  )

  list(
    A = A,
    Lambda = Lambda,
    X = X_res$X,
    Y = Y_res$X,
    c_r = c_r,
    c_c = c_c,
    Z_r = Z_r,
    Z_c = Z_c,
    M_r = X_res$M,
    M_c = Y_res$M,
    Sigma_r_list = X_res$Sigma_list,
    Sigma_c_list = Y_res$Sigma_list,
    model = "Poisson_BiDFM_no_degree_heterogeneity",
    rho = rho,
    B = B
  )
}

generate_connected_poisson_nodc_bisbm <- function(..., max_tries = 1000) {
  args <- list(...)
  for (attempt in seq_len(max_tries)) {
    dat <- do.call(generate_poisson_nodc_bisbm, args)
    if (all(rowSums(dat$A) > 0) && all(colSums(dat$A) > 0)) {
      dat$n_tries <- attempt
      return(dat)
    }
  }
  stop(sprintf(
    "在 %d 次尝试内未能生成无孤立节点的 Poisson 网络，请调大 rho 或 P 中的强度。",
    max_tries
  ))
}
