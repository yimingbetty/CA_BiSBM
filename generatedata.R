
pkgs <- c("MASS", "clue", "mclust", "ggplot2")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(MASS)
library(clue)
library(mclust)
library(ggplot2)

set.seed(2026)


OUTPUT_DIR <- "sim_outputs"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)


N_REPS_DEFAULT <- 10
T_GRID_DEFAULT <- 10


## =============================================================================
## 1. 协变量生成工具
## =============================================================================

#' 生成一个随机的正定协方差矩阵 (支持任意方向的各向异性)
#' @param p 维数
#' @param scale 特征值几何中心 (整体噪声水平)
#' @param condition 条件数 (最大/最小特征值), 1 表示各向同性, 越大越"扁"
generate_random_spd <- function(p, scale = 1, condition = 5) {
  if (p == 1) return(matrix(scale, 1, 1))
  Q <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  eigvals <- scale * exp(seq(-log(condition), log(condition), length.out = p))
  Sigma <- Q %*% diag(eigvals, p) %*% t(Q)
  (Sigma + t(Sigma)) / 2
}

#' 生成带社区结构 (可选组共享中心) 的协变量
#' @param cluster_groups 长度为 K 的母组编号向量; 同一母组下社区均值互相靠近
#'        (弱可分场景); 为 NULL 时每个社区独立采样中心
generate_covariates_with_groups <- function(n, K, labels, p,
                                            M = NULL,
                                            Sigma_list = NULL,
                                            sigma = 0.5,
                                            cluster_groups = NULL,
                                            mu_range = c(-2, 2),
                                            group_perturb_sd = 0.3) {
  if (is.null(M)) {
    if (!is.null(cluster_groups)) {
      n_groups <- max(cluster_groups)
      group_centers <- matrix(runif(n_groups * p, mu_range[1], mu_range[2]), n_groups, p)
      M <- matrix(0, K, p)
      for (k in 1:K) M[k, ] <- group_centers[cluster_groups[k], ] + rnorm(p, 0, group_perturb_sd)
    } else {
      M <- matrix(runif(K * p, mu_range[1], mu_range[2]), K, p)
    }
  }
  
  if (is.null(Sigma_list)) {
    Sigma_list <- lapply(1:K, function(k) diag(sigma^2, p))
  } else if (!is.list(Sigma_list)) {
    Sigma_single <- Sigma_list
    Sigma_list <- lapply(1:K, function(k) Sigma_single)
  }
  
  X <- matrix(NA_real_, n, p)
  for (i in seq_len(n)) X[i, ] <- MASS::mvrnorm(1, mu = M[labels[i], ], Sigma = Sigma_list[[labels[i]]])
  
  list(X = X, M = M, Sigma_list = Sigma_list)
}

#' 简单构造均值矩阵: 每个社区沿不同维度偏移 mu_gap (最简化的"可分度"控制方式)
#' 注意: 这个构造方式每个社区只在一个专属维度上有信号, 其余维度是纯噪声, 且
#' 社区 k 的信号强度正比于 k (社区1天生信号最弱)。仅作为快速演示用, 正式做
#' "协变量信号强度"相关的场景实验时建议改用下面的 make_full_signal_M()。
make_simple_M <- function(K, p, mu_gap) {
  M <- matrix(0, K, p)
  for (k in 1:K) M[k, ((k - 1) %% p) + 1] <- mu_gap * k
  M
}

#' 生成 K 个社区在 p 维空间中"稠密且均匀分离"的均值矩阵 (推荐: 比 make_simple_M 更合理)
#'
#' 与 make_simple_M 不同, 这里每个维度对所有社区都非零 (不再有"纯噪声维度"),
#' 且用信噪比 snr 直接、均匀地控制所有社区两两之间的可分度, 不会出现"某个社区
#' 天生信号更弱"的偏差。
#'
#' @param K 社区数
#' @param p 协变量维度
#' @param snr 目标信噪比: 任意两个社区中心之间的欧式距离 ≈ snr * sigma * sqrt(p)
#'        (sigma*sqrt(p) 近似是噪声在 p 维空间里的典型"半径")
#'        snr ≈ 0.5-1: 几乎不可分; snr ≈ 2-3: 中等可分度; snr >= 5: 很容易分开
#' @param sigma 协变量噪声标准差 (应与生成协变量时用的 sigma_r/sigma_c 保持一致)
make_full_signal_M <- function(K, p, snr = 2, sigma = 1) {
  raw <- matrix(rnorm(K * p), K, p)                 # 每个维度都非零, 不再稀疏
  raw <- scale(raw, center = TRUE, scale = FALSE)   # 让 K 个中心整体居中
  
  dists <- as.matrix(dist(raw))
  current_avg_dist <- mean(dists[upper.tri(dists)])
  target_dist <- snr * sigma * sqrt(p)
  scale_factor <- target_dist / current_avg_dist
  
  raw * scale_factor
}


## =============================================================================
## 2. Bernoulli CA-BiSBM 数据生成
## =============================================================================

#' 生成 CA-BiSBM 数据 (Bernoulli 网络 + 高斯协变量)
generate_ca_bisbm <- function(n_r, n_c, K_r, K_c, B,
                              p_r, p_c,
                              M_r = NULL, M_c = NULL,
                              Sigma_r_list = NULL, Sigma_c_list = NULL,
                              sigma_r = 1, sigma_c = 1,
                              row_cluster_groups = NULL, col_cluster_groups = NULL,
                              fixed = NULL) {
  
  if (!is.null(fixed) && !is.null(fixed$c_r) && !is.null(fixed$c_c)) {
    c_r <- fixed$c_r; c_c <- fixed$c_c
    Z_r <- if (!is.null(fixed$Z_r)) fixed$Z_r else { Zt <- matrix(0, n_r, K_r); Zt[cbind(1:n_r, c_r)] <- 1; Zt }
    Z_c <- if (!is.null(fixed$Z_c)) fixed$Z_c else { Zt <- matrix(0, n_c, K_c); Zt[cbind(1:n_c, c_c)] <- 1; Zt }
  } else {
    c_r <- sample(1:K_r, n_r, replace = TRUE)
    c_c <- sample(1:K_c, n_c, replace = TRUE)
    Z_r <- matrix(0, n_r, K_r); Z_r[cbind(1:n_r, c_r)] <- 1
    Z_c <- matrix(0, n_c, K_c); Z_c[cbind(1:n_c, c_c)] <- 1
  }
  
  Omega <- Z_r %*% B %*% t(Z_c)
  A <- matrix(rbinom(n_r * n_c, size = 1, prob = as.vector(Omega)), n_r, n_c)
  
  M_r_use     <- if (!is.null(fixed) && !is.null(fixed$M_r)) fixed$M_r else M_r
  M_c_use     <- if (!is.null(fixed) && !is.null(fixed$M_c)) fixed$M_c else M_c
  Sigma_r_use <- if (!is.null(fixed) && !is.null(fixed$Sigma_r_list)) fixed$Sigma_r_list else Sigma_r_list
  Sigma_c_use <- if (!is.null(fixed) && !is.null(fixed$Sigma_c_list)) fixed$Sigma_c_list else Sigma_c_list
  
  X_res <- generate_covariates_with_groups(n_r, K_r, c_r, p_r, M = M_r_use, Sigma_list = Sigma_r_use,
                                           sigma = sigma_r, cluster_groups = row_cluster_groups)
  Y_res <- generate_covariates_with_groups(n_c, K_c, c_c, p_c, M = M_c_use, Sigma_list = Sigma_c_use,
                                           sigma = sigma_c, cluster_groups = col_cluster_groups)
  
  list(A = A, X = X_res$X, Y = Y_res$X, c_r = c_r, c_c = c_c, Z_r = Z_r, Z_c = Z_c, Omega = Omega,
       M_r = X_res$M, M_c = Y_res$M, Sigma_r_list = X_res$Sigma_list, Sigma_c_list = Y_res$Sigma_list)
}

#' 保证无孤立节点的 CA-BiSBM 数据生成
generate_connected_ca_bisbm <- function(..., max_tries = 1000) {
  args <- list(...)
  for (attempt in 1:max_tries) {
    dat <- do.call(generate_ca_bisbm, args)
    if (all(rowSums(dat$A) > 0) && all(colSums(dat$A) > 0)) { dat$n_tries <- attempt; return(dat) }
  }
  stop(sprintf("在 %d 次尝试内未能生成无孤立节点的网络，请调大 B 中的边概率或减小 n_r, n_c。", max_tries))
}


## =============================================================================
## 3. 块概率矩阵 B 的构造 (支持 K_r != K_c)
## =============================================================================

#' 半自动构造 K_r x K_c 的块概率矩阵 (保证每个社区都至少有一个强连接对象)
make_B_asymmetric <- function(K_r, K_c, b_strong = 0.25, b_weak = 0.03) {
  B <- matrix(b_weak, K_r, K_c)
  if (K_r <= K_c) {
    assign_col <- rep(1:K_r, length.out = K_c)
    for (l in 1:K_c) B[assign_col[l], l] <- b_strong
  } else {
    assign_row <- rep(1:K_c, length.out = K_r)
    for (k in 1:K_r) B[k, assign_row[k]] <- b_strong
  }
  B
}


