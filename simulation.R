## =============================================================================
## 三种 B 矩阵下比较五种方法
##
## 1. Network only
## 2. Covariates only
## 3. NAC：Network-Adjusted Covariates
## 4. SDP：单模投影后的 Covariate-Regularized SDP
## 5. CA-BiSBM
##
## 网络模型：伯努利二部网络
## 本文件不修改 generatedata.R、methodnew.R 或之前的 simulator 文件。
## =============================================================================


source("generatedata.R")
source("generatedata_nonormal.R")
source("generatedata_poisson.R")
source("method.R")
library(ggplot2)

OUTPUT_DIR <- "sim_outputs"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

## -----------------------------------------------------------------------------
## 并行设置
## -----------------------------------------------------------------------------
HAS_FUTURE <- requireNamespace("future.apply", quietly = TRUE) &&
  requireNamespace("future", quietly = TRUE)

if (HAS_FUTURE) {
  library(future)
  library(future.apply)
  available_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(available_cores)) available_cores <- 2L
  N_WORKERS <- min(8L, max(1L, available_cores - 2L))
  future::plan(future::multisession, workers = N_WORKERS)
  cat(sprintf("并行 worker 数：%d\n", future::nbrOfWorkers()))
} else {
  message("future/future.apply 不可用，使用串行计算。")
}

## -----------------------------------------------------------------------------
## 仿真参数
## -----------------------------------------------------------------------------
N_REPS <- 100
T_GRID <- 20
NSTART <- 20

K_r <- 3
K_c <- 5
p_r <- 5
p_c <- 5
snr_r <- 1.5
snr_c <- 1.2
sigma_r <- 1
sigma_c <- 1

## 六个样本量点，保持 n_r:n_c = 3:5。
n_r_grid <- c(60, 120, 240, 480, 960, 1500)
n_c_grid <- round(n_r_grid * K_c / K_r)

## -----------------------------------------------------------------------------
## 三种伯努利 B 矩阵
## -----------------------------------------------------------------------------
B_balanced <- matrix(c(
  0.090, 0.030, 0.045, 0.065, 0.035,
  0.035, 0.090, 0.035, 0.040, 0.075,
  0.045, 0.035, 0.090, 0.060, 0.050
), nrow = 3, byrow = TRUE)

B_gradient <- matrix(c(
  0.100, 0.075, 0.045, 0.025, 0.015,
  0.020, 0.045, 0.100, 0.065, 0.025,
  0.015, 0.025, 0.045, 0.075, 0.100
), nrow = 3, byrow = TRUE)



B_TYPES <- list(
  balanced = B_balanced,
  gradient = B_gradient
)

if (any(vapply(B_TYPES, function(B) any(B < 0 | B > 1), logical(1)))) {
  stop("伯努利网络中的 B 元素必须全部位于 [0, 1] 内。")
}

## 固定协变量中心，使不同 B 场景之间只比较网络结构。
set.seed(2026)
M_r <- make_full_signal_M(K_r, p_r, snr = snr_r, sigma = sigma_r)
M_c <- make_full_signal_M(K_c, p_c, snr = snr_c, sigma = sigma_c)

## Case 4: Poisson intensity template without degree heterogeneity.
B_case4 <- matrix(c(
  2.1, 1.8, 1.9, 2.0, 1.8,
  1.9, 2.1, 1.8, 1.8, 2.0,
  1.8, 1.9, 2.1, 1.9, 1.9
), nrow = 3, byrow = TRUE)
RHO_CASE4 <- 0.8
n_r_grid_case4 <- c(60, 120, 240, 480, 960, 1500)
n_c_grid_case4 <- round(n_r_grid_case4 * K_c / K_r)

set.seed(2026)
M_r_case4 <- make_full_signal_M(K_r, p_r, snr = 1.8, sigma = sigma_r)
M_c_case4 <- make_full_signal_M(K_c, p_c, snr = 1.5, sigma = sigma_c)

## -----------------------------------------------------------------------------
## NAC：二部网络扩展
## -----------------------------------------------------------------------------
nac_one_side <- function(A_side, own_cov, neighbor_cov,
                         K_embed, K_target, nstart = 20) {
  n <- nrow(own_cov)
  degree <- rowSums(A_side)
  log_n <- log(max(n, 2))
  
  alpha <- (mean(degree) / 2) * log_n / (degree + log_n)
  adjusted_cov <- A_side %*% neighbor_cov + own_cov * alpha
  
  K_use <- min(K_embed, ncol(adjusted_cov), n - 1L)
  sv <- svd(adjusted_cov, nu = K_use, nv = 0)
  U <- sv$u[, seq_len(K_use), drop = FALSE]
  U <- U / pmax(sqrt(rowSums(U^2)), .Machine$double.eps)
  
  kmeans(U, centers = K_target, nstart = nstart, iter.max = 100)$cluster
}

nac_bipartite <- function(A, X, Y, K_r, K_c, nstart = 20) {
  K_embed <- min(K_r, K_c)
  list(
    c_r_hat = nac_one_side(A, X, Y, K_embed, K_r, nstart),
    c_c_hat = nac_one_side(t(A), Y, X, K_embed, K_c, nstart)
  )
}

## -----------------------------------------------------------------------------
## SDP：先投影为两个单模网络，再使用 NAC::SDP
## -----------------------------------------------------------------------------
SDP_LAMBDA <- 0.20
SDP_ALPHA <- 0.50
SDP_RHO <- 2
SDP_TT <- 100
SDP_TOL <- 5

make_projection_adjacency <- function(A, side = c("row", "column")) {
  side <- match.arg(side)
  S <- if (side == "row") A %*% t(A) else t(A) %*% A
  Adj <- 1L * (S > 0)
  diag(Adj) <- 0L
  Adj
}

run_sdp <- function(Adj, Covariate, K) {
  result <- NULL
  invisible(capture.output({
    result <- NAC::SDP(
      Adj = Adj,
      Covariate = Covariate,
      lambda = SDP_LAMBDA,
      K = K,
      alpha = SDP_ALPHA,
      rho = SDP_RHO,
      TT = SDP_TT,
      tol = SDP_TOL,
      quiet = TRUE
    )
  }))
  as.integer(result)
}

sdp_bipartite_projection <- function(A, X, Y, K_r, K_c) {
  row_labels <- run_sdp(
    make_projection_adjacency(A, "row"), X, K_r
  )
  col_labels <- run_sdp(
    make_projection_adjacency(A, "column"), Y, K_c
  )
  list(c_r_hat = row_labels, c_c_hat = col_labels)
}

## -----------------------------------------------------------------------------
## 单次实验
## -----------------------------------------------------------------------------
generate_one_dat <- function(B, n_r, n_c, seed,
                             case_name = "bernoulli") {
  set.seed(seed)
  
  if (case_name == "nonnormal") {
    return(generate_connected_nonnormal_ca_bisbm(
      n_r = n_r, n_c = n_c, K_r = K_r, K_c = K_c,
      B = B, p_r = p_r, p_c = p_c,
      M_r = M_r, M_c = M_c,
      sigma_r = sigma_r, sigma_c = sigma_c,
      df_r = 3, df_c = 3, max_tries = 1000
    ))
  }
  
  if (case_name == "poisson") {
    return(generate_connected_poisson_nodc_bisbm(
      n_r = n_r, n_c = n_c, K_r = K_r, K_c = K_c,
      B = B, rho = RHO_CASE4,
      p_r = p_r, p_c = p_c,
      M_r = M_r_case4, M_c = M_c_case4,
      sigma_r = sigma_r, sigma_c = sigma_c,
      max_tries = 1000
    ))
  }
  
  generate_connected_ca_bisbm(
    n_r = n_r, n_c = n_c, K_r = K_r, K_c = K_c,
    B = B, p_r = p_r, p_c = p_c,
    M_r = M_r, M_c = M_c,
    sigma_r = sigma_r, sigma_c = sigma_c,
    max_tries = 1000
  )
}

run_one_main <- function(B, n_r, n_c, seed, case_name = "bernoulli") {
  dat <- generate_one_dat(B, n_r, n_c, seed, case_name)
  
  base_fit <- run_methods(
    dat = dat,
    K_r = K_r,
    K_c = K_c,
    T_grid = T_GRID,
    nstart = NSTART
  )
  
  nac_fit <- nac_bipartite(
    A = dat$A,
    X = dat$X,
    Y = dat$Y,
    K_r = K_r,
    K_c = K_c,
    nstart = NSTART
  )
  
  data.frame(
    err_net_r = base_fit$err_net_r,
    err_cov_r = base_fit$err_cov_r,
    err_ca_r = base_fit$err_ca_r,
    err_nac_r = clustering_error(dat$c_r, nac_fit$c_r_hat),
    err_net_c = base_fit$err_net_c,
    err_cov_c = base_fit$err_cov_c,
    err_ca_c = base_fit$err_ca_c,
    err_nac_c = clustering_error(dat$c_c, nac_fit$c_c_hat),
    ari_net_r = base_fit$ari_net_r,
    ari_cov_r = base_fit$ari_cov_r,
    ari_ca_r = base_fit$ari_ca_r,
    ari_nac_r = mclust::adjustedRandIndex(dat$c_r, nac_fit$c_r_hat),
    ari_net_c = base_fit$ari_net_c,
    ari_cov_c = base_fit$ari_cov_c,
    ari_ca_c = base_fit$ari_ca_c,
    ari_nac_c = mclust::adjustedRandIndex(dat$c_c, nac_fit$c_c_hat),
    alpha_hat = base_fit$alpha_hat,
    beta_hat = base_fit$beta_hat
  )
}

run_one_sdp <- function(B, n_r, n_c, seed, case_name = "bernoulli") {
  dat <- generate_one_dat(B, n_r, n_c, seed, case_name)
  
  sdp_fit <- sdp_bipartite_projection(
    A = dat$A,
    X = dat$X,
    Y = dat$Y,
    K_r = K_r,
    K_c = K_c
  )
  
  data.frame(
    err_sdp_r = clustering_error(dat$c_r, sdp_fit$c_r_hat),
    err_sdp_c = clustering_error(dat$c_c, sdp_fit$c_c_hat),
    ari_sdp_r = mclust::adjustedRandIndex(dat$c_r, sdp_fit$c_r_hat),
    ari_sdp_c = mclust::adjustedRandIndex(dat$c_c, sdp_fit$c_c_hat)
  )
}

run_reps <- function(B, n_r, n_c, n_reps, worker_fun,
                     case_name = "bernoulli", seed_base = 9000) {
  worker <- function(r) {
    worker_fun(B, n_r, n_c, seed_base + r, case_name)
  }
  
  if (HAS_FUTURE && n_reps > 1L) {
    reps <- tryCatch(
      future.apply::future_lapply(
        seq_len(n_reps), worker, future.seed = TRUE
      ),
      error = function(e) {
        warning("并行计算失败，切换为串行：", conditionMessage(e))
        future::plan(future::sequential)
        lapply(seq_len(n_reps), worker)
      }
    )
  } else {
    reps <- lapply(seq_len(n_reps), worker)
  }
  
  do.call(rbind, reps)
}

make_replicate_long <- function(main_reps, sdp_reps, B_name, n_r, n_c) {
  n_reps <- nrow(main_reps)
  rep_id <- seq_len(n_reps)
  
  rbind(
    data.frame(
      B_type = B_name, n_r = n_r, n_c = n_c, rep = rep_id,
      method = "Network only",
      row_error = main_reps$err_net_r,
      column_error = main_reps$err_net_c,
      row_ARI = main_reps$ari_net_r,
      column_ARI = main_reps$ari_net_c,
      alpha_hat = NA_real_, beta_hat = NA_real_
    ),
    data.frame(
      B_type = B_name, n_r = n_r, n_c = n_c, rep = rep_id,
      method = "Covariates only",
      row_error = main_reps$err_cov_r,
      column_error = main_reps$err_cov_c,
      row_ARI = main_reps$ari_cov_r,
      column_ARI = main_reps$ari_cov_c,
      alpha_hat = NA_real_, beta_hat = NA_real_
    ),
    data.frame(
      B_type = B_name, n_r = n_r, n_c = n_c, rep = rep_id,
      method = "NAC",
      row_error = main_reps$err_nac_r,
      column_error = main_reps$err_nac_c,
      row_ARI = main_reps$ari_nac_r,
      column_ARI = main_reps$ari_nac_c,
      alpha_hat = NA_real_, beta_hat = NA_real_
    ),
    data.frame(
      B_type = B_name, n_r = n_r, n_c = n_c, rep = rep_id,
      method = "CA-BiSBM",
      row_error = main_reps$err_ca_r,
      column_error = main_reps$err_ca_c,
      row_ARI = main_reps$ari_ca_r,
      column_ARI = main_reps$ari_ca_c,
      alpha_hat = main_reps$alpha_hat,
      beta_hat = main_reps$beta_hat
    ),
    data.frame(
      B_type = B_name, n_r = n_r, n_c = n_c, rep = seq_len(nrow(sdp_reps)),
      method = "SDP",
      row_error = sdp_reps$err_sdp_r,
      column_error = sdp_reps$err_sdp_c,
      row_ARI = sdp_reps$ari_sdp_r,
      column_ARI = sdp_reps$ari_sdp_c,
      alpha_hat = NA_real_, beta_hat = NA_real_
    )
  )
}

## -----------------------------------------------------------------------------
## 运行三个 B 场景
## -----------------------------------------------------------------------------
cat("\n########## 三种 B 矩阵下的五方法比较 ##########\n")
print(data.frame(n_r = n_r_grid, n_c = n_c_grid))

replicate_results <- list()

CASE_CONFIGS <- list(
  scenario1_distinct = list(B = B_gradient, case_name = "bernoulli"),
  scenario2_balanced = list(B = B_balanced, case_name = "bernoulli"),
  scenario3_nonnormal = list(B = B_balanced, case_name = "nonnormal"),
  scenario4_poisson = list(B = B_case4, case_name = "poisson")
)

comparison_results <- do.call(rbind, lapply(names(CASE_CONFIGS), function(B_name) {
  config <- CASE_CONFIGS[[B_name]]
  B <- config$B
  case_name <- config$case_name
  
  do.call(rbind, Map(function(nr, nc) {
    cat(sprintf("Scenario = %s, n_r = %d, n_c = %d\n", B_name, nr, nc))
    
    main_reps <- run_reps(
      B = B,
      n_r = nr,
      n_c = nc,
      n_reps = N_REPS,
      worker_fun = run_one_main,
      case_name = case_name,
      seed_base = 9000 + nr
    )
    
    sdp_reps <- run_reps(
      B = B,
      n_r = nr,
      n_c = nc,
      n_reps = N_REPS,
      worker_fun = run_one_sdp,
      case_name = case_name,
      seed_base = if (case_name == "poisson") 9000 + nr else 19000 + nr
    )
    
    replicate_results[[length(replicate_results) + 1L]] <<-
      make_replicate_long(
        main_reps = main_reps,
        sdp_reps = sdp_reps,
        B_name = B_name,
        n_r = nr,
        n_c = nc
      )
    
    data.frame(
      B_type = B_name,
      n_r = nr,
      n_c = nc,
      err_net_r = mean(main_reps$err_net_r),
      err_cov_r = mean(main_reps$err_cov_r),
      err_ca_r = mean(main_reps$err_ca_r),
      err_nac_r = mean(main_reps$err_nac_r),
      err_sdp_r = mean(sdp_reps$err_sdp_r),
      err_net_c = mean(main_reps$err_net_c),
      err_cov_c = mean(main_reps$err_cov_c),
      err_ca_c = mean(main_reps$err_ca_c),
      err_nac_c = mean(main_reps$err_nac_c),
      err_sdp_c = mean(sdp_reps$err_sdp_c),
      ari_net_r = mean(main_reps$ari_net_r),
      ari_cov_r = mean(main_reps$ari_cov_r),
      ari_ca_r = mean(main_reps$ari_ca_r),
      ari_nac_r = mean(main_reps$ari_nac_r),
      ari_sdp_r = mean(sdp_reps$ari_sdp_r),
      ari_net_c = mean(main_reps$ari_net_c),
      ari_cov_c = mean(main_reps$ari_cov_c),
      ari_ca_c = mean(main_reps$ari_ca_c),
      ari_nac_c = mean(main_reps$ari_nac_c),
      ari_sdp_c = mean(sdp_reps$ari_sdp_c),
      alpha_hat = mean(main_reps$alpha_hat),
      beta_hat = mean(main_reps$beta_hat)
    )
  }, if (case_name == "poisson") n_r_grid_case4 else n_r_grid,
  if (case_name == "poisson") n_c_grid_case4 else n_c_grid))
}))

print(comparison_results)
replicate_results_df <- do.call(rbind, replicate_results)
write.csv(
  comparison_results,
  file.path(OUTPUT_DIR, "comparison_four_cases_five_methods_ARI.csv"),
  row.names = FALSE
)
write.csv(
  replicate_results_df,
  file.path(OUTPUT_DIR, "comparison_four_cases_five_methods_ARI_replicates.csv"),
  row.names = FALSE
)

## -----------------------------------------------------------------------------
## 绘图
## -----------------------------------------------------------------------------
plot_df <- rbind(
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "Network only", error = comparison_results$err_net_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "Covariates only", error = comparison_results$err_cov_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "NAC", error = comparison_results$err_nac_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "SDP", error = comparison_results$err_sdp_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "CA-BiSBM", error = comparison_results$err_ca_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "Network only", error = comparison_results$err_net_c),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "Covariates only", error = comparison_results$err_cov_c),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "NAC", error = comparison_results$err_nac_c),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "SDP", error = comparison_results$err_sdp_c),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "CA-BiSBM", error = comparison_results$err_ca_c)
)

plot_df$method <- factor(
  plot_df$method,
  levels = c("Network only", "Covariates only", "NAC", "SDP", "CA-BiSBM")
)

p <- ggplot(plot_df, aes(x = n_r, y = error, color = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.7) +
  facet_grid(B_type ~ side) +
  scale_x_continuous(breaks = n_r_grid) +
  labs(
    x = expression(n[r]),
    y = "Average misclassification rate",
    title = "Five-method comparison under four simulation scenarios",
    subtitle = "Bernoulli, non-Gaussian, and Poisson bipartite networks"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  file.path(OUTPUT_DIR, "comparison_four_cases_five_methods_ARI_errors.png"),
  p, width = 12, height = 8, dpi = 150
)

## -----------------------------------------------------------------------------
## ARI plot
## -----------------------------------------------------------------------------
ari_df <- rbind(
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "Network only", ARI = comparison_results$ari_net_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "Covariates only", ARI = comparison_results$ari_cov_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "NAC", ARI = comparison_results$ari_nac_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "SDP", ARI = comparison_results$ari_sdp_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Row nodes", method = "CA-BiSBM", ARI = comparison_results$ari_ca_r),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "Network only", ARI = comparison_results$ari_net_c),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "Covariates only", ARI = comparison_results$ari_cov_c),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "NAC", ARI = comparison_results$ari_nac_c),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "SDP", ARI = comparison_results$ari_sdp_c),
  data.frame(B_type = comparison_results$B_type, n_r = comparison_results$n_r, side = "Column nodes", method = "CA-BiSBM", ARI = comparison_results$ari_ca_c)
)

ari_df$method <- factor(
  ari_df$method,
  levels = c("Network only", "Covariates only", "NAC", "SDP", "CA-BiSBM")
)

p_ari <- ggplot(ari_df, aes(x = n_r, y = ARI, color = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.7) +
  facet_grid(B_type ~ side) +
  scale_x_continuous(breaks = n_r_grid) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = expression(n[r]),
    y = "Average adjusted Rand index",
    title = "ARI comparison under four simulation scenarios",
    subtitle = "Bernoulli, non-Gaussian, and Poisson bipartite networks"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  file.path(OUTPUT_DIR, "comparison_four_cases_five_methods_ARI.png"),
  p_ari, width = 12, height = 8, dpi = 150
)

cat("\n已保存：sim_outputs/comparison_three_B_five_methods.csv\n")
