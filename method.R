

if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")
if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
library(clue)
library(mclust)



construct_L <- function(A) {
  d_r <- rowSums(A); d_c <- colSums(A)
  d_bar_r <- mean(d_r); d_bar_c <- mean(d_c)
  D_r_inv_sqrt <- 1 / sqrt(d_r + d_bar_r)
  D_c_inv_sqrt <- 1 / sqrt(d_c + d_bar_c)
  L <- sweep(A, 1, D_r_inv_sqrt, "*")
  L <- sweep(L, 2, D_c_inv_sqrt, "*")
  L
}

normalize_rows <- function(U) {
  norms <- sqrt(rowSums(U^2)); norms[norms == 0] <- 1e-10
  U / norms
}

center_covariates <- function(X) {
  scale(X, center = TRUE, scale = FALSE)
}

normalize_psd <- function(S) {
  top_eigenvalue <- max(
    eigen(S, symmetric = TRUE, only.values = TRUE)$values
  )
  if (!is.finite(top_eigenvalue) || top_eigenvalue <= 1e-12) {
    return(S)
  }
  S / top_eigenvalue
}

make_log_grid <- function(lo, hi, n) {
  lo <- max(abs(lo), 1e-8)
  hi <- max(abs(hi), lo * 1.01)
  exp(seq(log(lo), log(hi), length.out = n))
}

ca_spectral_clustering <- function(L, X, Y, K_r, K_c, alpha, beta, nstart = 20) {
  K <- min(K_r, K_c)

  network_r <- L %*% t(L)
  network_c <- t(L) %*% L
  cov_r <- X %*% t(X)
  cov_c <- Y %*% t(Y)

  L_r <- network_r + alpha * cov_r
  L_c <- network_c + beta * cov_c
  
  eig_r <- eigen(L_r, symmetric = TRUE)
  eig_c <- eigen(L_c, symmetric = TRUE)
  
  U_r_star <- normalize_rows(eig_r$vectors[, 1:K, drop = FALSE])
  U_c_star <- normalize_rows(eig_c$vectors[, 1:K, drop = FALSE])
  
  km_r <- kmeans(U_r_star, centers = K_r, nstart = nstart, iter.max = 100)
  km_c <- kmeans(U_c_star, centers = K_c, nstart = nstart, iter.max = 100)
  
  list(c_r_hat = km_r$cluster, c_c_hat = km_c$cluster,
       J_r = km_r$tot.withinss, J_c = km_c$tot.withinss)
}

select_tuning_params <- function(L, X, Y, K_r, K_c, T_grid = 10, nstart = 10) {
  K <- min(K_r, K_c)
  L_r_base <- L %*% t(L)
  L_c_base <- t(L) %*% L
  XXt <- X %*% t(X)
  YYt <- Y %*% t(Y)
  if (length(T_grid) != 1L || !is.finite(T_grid) || T_grid < 2L || T_grid != as.integer(T_grid)) {
    stop("T_grid must be an integer >= 2 for the paper's equally spaced grid.")
  }
  if (K < 1L || K >= min(nrow(L), ncol(L))) {
    stop("K must allow both the Kth and (K+1)th eigenvalues.")
  }

  eig_LLt <- eigen(L_r_base, symmetric = TRUE, only.values = TRUE)$values
  eig_LtL <- eigen(L_c_base, symmetric = TRUE, only.values = TRUE)$values
  eig_XXt <- eigen(XXt, symmetric = TRUE, only.values = TRUE)$values
  eig_YYt <- eigen(YYt, symmetric = TRUE, only.values = TRUE)$values
  p_r <- ncol(X); p_c <- ncol(Y)
  
  alpha_min <- (eig_LLt[K] - eig_LLt[K + 1]) / eig_XXt[1]
  denom_a <- if (p_r <= K) eig_XXt[p_r] else (eig_XXt[K] - eig_XXt[K + 1])
  alpha_max <- eig_LLt[1] / denom_a
  
  beta_min <- (eig_LtL[K] - eig_LtL[K + 1]) / eig_YYt[1]
  denom_b <- if (p_c <= K) eig_YYt[p_c] else (eig_YYt[K] - eig_YYt[K + 1])
  beta_max <- eig_LtL[1] / denom_b
  
  # Invalid/degenerate intervals are not repaired by changing the paper formula.
  fix_interval <- function(lo, hi) {
    if (any(!is.finite(c(lo, hi))) || lo < 0 || hi < lo) {
      stop("Paper tuning interval is invalid; inspect covariate rank and eigengaps.")
    }
    c(lo, hi)
  }
  ab_a <- fix_interval(alpha_min, alpha_max)
  ab_b <- fix_interval(beta_min, beta_max)
  


  
  alpha_grid <- c(0, make_log_grid(ab_a[1], ab_a[2], T_grid))
  beta_grid  <- c(0, make_log_grid(ab_b[1], ab_b[2], T_grid))
  
  row_fits <- lapply(alpha_grid, function(a) {
    eig_r <- eigen(L_r_base + a * XXt, symmetric = TRUE)
    U <- normalize_rows(eig_r$vectors[, 1:K, drop = FALSE])
    kmeans(U, centers = K_r, nstart = nstart, iter.max = 100)
  })
  column_fits <- lapply(beta_grid, function(b) {
    eig_c <- eigen(L_c_base + b * YYt, symmetric = TRUE)
    U <- normalize_rows(eig_c$vectors[, 1:K, drop = FALSE])
    kmeans(U, centers = K_c, nstart = nstart, iter.max = 100)
  })
  
  J_r_vec <- vapply(row_fits, function(z) z$tot.withinss, numeric(1))
  J_c_vec <- vapply(column_fits, function(z) z$tot.withinss, numeric(1))
  ir <- which.min(J_r_vec); ic <- which.min(J_c_vec)
  list(alpha_hat = alpha_grid[ir], beta_hat = beta_grid[ic],
       alpha_range = ab_a, beta_range = ab_b,
       alpha_grid = alpha_grid, beta_grid = beta_grid,
       J_r_vec = J_r_vec, J_c_vec = J_c_vec,
       fit = list(c_r_hat = row_fits[[ir]]$cluster,
                  c_c_hat = column_fits[[ic]]$cluster,
                  J_r = J_r_vec[ir], J_c = J_c_vec[ic]))
}



select_tuning_params_oracle <- function(
    L, X, Y, c_r_true, c_c_true,
    K_r, K_c, T_grid = 10
) {
  K <- min(K_r, K_c)
  
  LLt <- L %*% t(L)
  LtL <- t(L) %*% L
  XXt <- X %*% t(X)
  YYt <- Y %*% t(Y)
  
  tmp <- select_tuning_params(
    L, X, Y,
    K_r = K_r,
    K_c = K_c,
    T_grid = T_grid
  )
  

  alpha_grid <- tmp$alpha_grid
  beta_grid <- tmp$beta_grid
  
  get_embedding <- function(S) {
    eig <- eigen(S, symmetric = TRUE)
    
    U <- eig$vectors[, 1:K, drop = FALSE]
    
    normalize_rows(U)
  }
  
  err_r_vec <- sapply(alpha_grid, function(a) {
    U_r <- get_embedding(LLt + a * XXt)
    
    km_r <- kmeans(
      U_r,
      centers = K_r,
      nstart = 30,
      iter.max = 100
    )
    
    clustering_error(c_r_true, km_r$cluster)
  })
  
  err_c_vec <- sapply(beta_grid, function(b) {
    U_c <- get_embedding(LtL + b * YYt)
    
    km_c <- kmeans(
      U_c,
      centers = K_c,
      nstart = 30,
      iter.max = 100
    )
    
    clustering_error(c_c_true, km_c$cluster)
  })
  
  list(
    alpha_hat = alpha_grid[which.min(err_r_vec)],
    beta_hat = beta_grid[which.min(err_c_vec)],
    best_err_r = min(err_r_vec),
    best_err_c = min(err_c_vec),
    alpha_grid = alpha_grid,
    beta_grid = beta_grid,
    err_r_vec = err_r_vec,
    err_c_vec = err_c_vec
  )
}



clustering_error <- function(true_labels, est_labels) {
  K <- max(length(unique(true_labels)), length(unique(est_labels)))
  conf_mat <- matrix(0, K, K)
  for (k in 1:K) for (l in 1:K) conf_mat[k, l] <- sum(true_labels == k & est_labels == l)
  assignment <- clue::solve_LSAP(conf_mat, maximum = TRUE)
  1 - sum(conf_mat[cbind(1:K, assignment)]) / length(true_labels)
}

## Covariate-only baseline using the same K=min(K_r,K_c) embedding dimension
## as CA-BiSBM, for a fair comparison when K_r != K_c.
covariate_spectral_clustering <- function(W, K_embed, K_target, nstart = 20) {
  S <- W %*% t(W)

  eig <- eigen(S, symmetric = TRUE)
  U <- eig$vectors[, 1:K_embed, drop = FALSE]
  U <- normalize_rows(U)

  kmeans(
    U,
    centers = K_target,
    nstart = nstart,
    iter.max = 100
  )
}



run_methods <- function(dat, K_r, K_c, T_grid = T_GRID_DEFAULT, nstart = 15) {
  
  L <- construct_L(dat$A)
  
  # (a) 仅网络: alpha = beta = 0
  res_net <- ca_spectral_clustering(L, dat$X, dat$Y, K_r, K_c, alpha = 0, beta = 0, nstart = nstart)
  
  # (b) 仅协变量: 对原始 Gram 矩阵的行归一化谱嵌入做 K-means
  K <- min(K_r, K_c)
  km_cov_r <- covariate_spectral_clustering(
    dat$X, K_embed = K, K_target = K_r, nstart = nstart
  )
  km_cov_c <- covariate_spectral_clustering(
    dat$Y, K_embed = K, K_target = K_c, nstart = nstart
  )
  
  # (c) CA-BiSBM: Algorithm 2 数据驱动选取 alpha, beta, 再做 Algorithm 1
  tune <- select_tuning_params(L, dat$X, dat$Y, K_r, K_c, T_grid = T_grid, nstart = min(nstart, 10))
  # Paper p. 7: retain partitions from the winning WCSS candidates.
  res_ca <- tune$fit
  
  data.frame(
    err_net_r = clustering_error(dat$c_r, res_net$c_r_hat),
    err_net_c = clustering_error(dat$c_c, res_net$c_c_hat),
    err_cov_r = clustering_error(dat$c_r, km_cov_r$cluster),
    err_cov_c = clustering_error(dat$c_c, km_cov_c$cluster),
    err_ca_r  = clustering_error(dat$c_r, res_ca$c_r_hat),
    err_ca_c  = clustering_error(dat$c_c, res_ca$c_c_hat),
    ari_net_r = mclust::adjustedRandIndex(dat$c_r, res_net$c_r_hat),
    ari_net_c = mclust::adjustedRandIndex(dat$c_c, res_net$c_c_hat),
    ari_cov_r = mclust::adjustedRandIndex(dat$c_r, km_cov_r$cluster),
    ari_cov_c = mclust::adjustedRandIndex(dat$c_c, km_cov_c$cluster),
    ari_ca_r  = mclust::adjustedRandIndex(dat$c_r, res_ca$c_r_hat),
    ari_ca_c  = mclust::adjustedRandIndex(dat$c_c, res_ca$c_c_hat),
    alpha_hat = tune$alpha_hat, beta_hat = tune$beta_hat
  )
}
