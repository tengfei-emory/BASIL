
compute_jic <- function(Y, svd_Y, k){
  n <- nrow(Y); p <- ncol(Y)
  minint <- min(n ,p)
  maxint <- max(n, p)
  M <- sqrt(n)*as.matrix(svd_Y$u[,1:k])
  Lambda <- 1/sqrt(n)* as.matrix(svd_Y$v[,1:k]) %*% diag(svd_Y$d[1:k], k, k)
  Y_hat <- tcrossprod(M, Lambda)
  sigma_sq_hat <- colMeans((Y - Y_hat)^2) # p * 1
  tausq_est <- (mean(colSums((Y_hat)^2) / sigma_sq_hat)) / (n * k);
  Lambda_est <- (sqrt(n) / (n + 1/tausq_est)) * svd_Y$v[,1:k] %*% diag(svd_Y$d[1:k], k, k);
  M_Lambda_est <- sqrt(n)* svd_Y$u[,1:k] %*% t(Lambda_est)
  res <- colSums((Y - M_Lambda_est)^2)
  sigma_sq_hat <- res / n
  term1 <- -n*sum(log(sqrt(sigma_sq_hat))) - 0.5*sum(res/sigma_sq_hat)
  term1 <- term1 / (n*p)
  jic <- -2*term1 + (k * maxint * log(minint) / (n*p))

  return(jic)
}

compute_posterior_samples <- function(
    Y, Lambda_C, Lambda_N, tau_C, tau_N, sigma_sq, U_C,
    v0=1, sigma_sq_0=1, n_MC=100
){

  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(Lambda_C)
  q <- ncol(U_C)

  SST <- sum(Y^2)
  SSE <- (n + 1/tau_C) * sum(Lambda_C^2) + (n + 1/tau_N) * sum(Lambda_N^2)
  sigma_sq_save <- rep(NA, n_MC)
  v_n <- v0 + n*p
  gamma_n <- v0*sigma_sq_0 + SST - SSE
  sigma_sq_save <- 1/rgamma(n_MC, v_n/2, gamma_n/2)

  Lambda_save <- array(NA, dim=c(p, k, n_MC))

  if(tau_N>tau_C){
    P_C <- tcrossprod(U_C)
    Q_C <- diag(1, p, p)- P_C
    N_C <- svd(Q_C)$u[,1:(p-q)]
  }

  for(t in 1:n_MC){
    if(tau_N<tau_C){
      E <- matrix(rnorm(p*k), ncol=k) * sqrt(sigma_sq_save[t] / (n + 1/tau_N))
      E_2 <- U_C %*% matrix(rnorm(q*k), ncol=k) * sqrt(sigma_sq_save[t] * (1/tau_N - 1/tau_C) / (n + 1/tau_C))
    } else {
      E <- matrix(rnorm(p*k), ncol=k) * sqrt(sigma_sq_save[t] / (n + 1/tau_C))
      E_2 <- N_C %*% matrix(rnorm((p-q)*k), ncol=k) * sqrt(sigma_sq_save[t] * (1/tau_C - 1/tau_N) / (n + 1/tau_N))
    }
    Lambda_save[,,t] <- Lambda_C + Lambda_N + E + E_2
  }
  return(list(Lambda_samples = Lambda_save, sigma_sq_samples = sigma_sq_save))
}
