
#' Estimate the latent dimension using Joint Information Criterion
#' @param Y A numeric matrix of gene expression data (samples × genes)
#' @param k_max Maximum number of latent factors to consider.
#'
#' @return A list containing:
#'   \item{k_hat}{Integer, the estimated optimal number of latent factors}
#'   \item{jics}{Numeric vector of JIC values for k = 1 to k_max}
#'   \item{svd_Y}{SVD decomposition of Y (output from svd() function)}
#'
#' @export
estimate_latent_dimension <- function(Y, k_max){
  svd_Y <- svd(Y)
  n <- nrow(Y)
  jics <- sapply(1:k_max, function(x) (compute_jic(Y, svd_Y, x)))
  plot(1:k_max, jics, type='l', xlab='k', ylab='jic', main='')
  print(paste('k_hat = ', which.min(jics)))
  return(list(k_hat = which.min(jics), jics=jics, svd_Y = svd_Y))
}


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

#' Compute point estimates for BASIL model parameters
#'
#' This function computes point estimates of the loading matrices (Lambda_C and
#' Lambda_N), latent factors (M), error variance (sigma_sq), and prior precision
#' parameters (tau_C and tau_N) for the pathway-based Bayesian factor model.
#'
#' @param Y A numeric matrix of gene expression data (samples × genes)
#' @param C A binary gene membership matrix (genes × pathways).
#' @param k Integer specifying the number of latent factors. If NA (default),
#'   the function will estimate k using \code{estimate_latent_dimension}.
#' @param k_max Maximum number of factors to consider when estimating k.
#'   Only used if k = NA. Default is 50.
#' @param v0 Prior degrees of freedom for the error variance. Default is 1.
#' @param sigma_sq0 Prior scale parameter for the error variance. Default is 1.
#'
#' @return A list containing:
#'   \item{M}{Matrix of latent factor scores (samples × k)}
#'   \item{Lambda_C}{Loading matrix for common/pathway-guided factors (genes × k)}
#'   \item{Lambda_N}{Loading matrix for noise/non-pathway factors (genes × k)}
#'   \item{sigma_sq}{Estimated error variance (scalar)}
#'   \item{tau_C}{Prior precision parameter for common factors (scalar)}
#'   \item{tau_N}{Prior precision parameter for noise factors (scalar)}
#'   \item{k}{Number of latent factors used (integer)}
#'   \item{P_C}{Projection matrix onto pathway space (genes × genes)}
#'
#' @export
#'
compute_point_estimates <- function(
    Y, C, k = NA, k_max = 50, v0 = 1, sigma_sq0 = 1){

  n <- nrow(Y); p <- ncol(Y)

  # estimate latent dimension if unknown
  if(is.na(k)){
    k_hat <-  estimate_latent_dimension(Y, k_max)
    k <- k_hat$k_hat
    svd_Y <- k_hat$svd_Y
  }
  else{
    svd_Y <- svd(Y)
  }

  # extract components of svd of Y
  U <- svd_Y$u[,1:k]
  D <- diag(as.vector(svd_Y$d[1:k]))
  V <- svd_Y$v[,1:k]
  U_perp <- svd_Y$u[,-c(1:k)]
  D_perp <- diag(as.vector(svd_Y$d[-c(1:k)]))
  V_perp <- svd_Y$v[,-c(1:k)]

  #P_C <- C %*% solve(crossprod(C)) %*% t(C)
  tCCt <- crossprod(C)
  s_cross_C <- svd(tCCt)
  q <- sum(s_cross_C$d>(0.1)^7)
  V_C <- s_cross_C$u[,1:q]
  D_C <- diag(as.vector(sqrt(s_cross_C$d[1:q])))
  D_C_inv <- diag(as.vector(1/sqrt(s_cross_C$d[1:q])))
  U_C <- C %*% V_C %*% D_C_inv
  P_C <- tcrossprod(U_C)
  Q_C <- diag(1, p, p) - P_C

  Y_P <- Y %*% P_C
  Y_Q <- Y - Y_P
  V_D <- V %*% D
  P_V_D <- P_C %*% V_D
  Q_V_D <- V_D - P_V_D
  D_Vt_perp <- D_perp %*% t(V_perp)
  D_Vt_perp_P <-  D_Vt_perp %*% P_C


  sigma_sq_hat <- sum(D_perp^2) / ((n - k)*p)

  # prior variances hyperparms
  #tau_C <- sum((P_V_D)^2 ) / (k * sum((D_Vt_perp_P )^2) ) * (n-k) / n
  tau_C <- sum((P_V_D)^2 ) / (k * sigma_sq_hat * q) / n
  print(paste('tau_C  = ', tau_C))
  #tau_N <- sum((Q_V_D)^2 ) / (k * sum((D_Vt_perp - D_Vt_perp_P)^2) ) * (n-k) / n
  tau_N <- sum((Q_V_D)^2 ) / (k * sigma_sq_hat * (p-q))/ n
  print(paste('tau_N  = ', tau_N))

  # point estimates
  M <- sqrt(n) * U
  Lambda_C <- (P_V_D) * sqrt(n) / (n + 1/tau_C)
  Lambda_N <- (Q_V_D) * sqrt(n) / (n + 1/tau_N)
  sigma_sq <- (v0 * sigma_sq0 + sum((Y-tcrossprod(M, Lambda_C + Lambda_N) )^2) ) /
    (v0 + n*p -2)

  return(list(M=M, Lambda_C=Lambda_C, Lambda_N=Lambda_N, sigma_sq=sigma_sq,
              tau_C=tau_C, tau_N=tau_N, k=k, P_C=P_C))
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



#' Compute posterior samples with correlation correction
#'
#' @param Y A numeric matrix of gene expression data (samples × genes)
#' @param Lambda_C Loading matrix for common/pathway-guided factors (genes × k)
#' @param Lambda_N Loading matrix for noise/non-pathway factors (genes × k)
#' @param tau_C Prior precision parameter for common factors
#' @param tau_N Prior precision parameter for noise factors
#' @param sigma_sq Error variance estimate
#' @param P_C Projection matrix onto pathway space (genes × genes)
#' @param v0 Prior degrees of freedom for error variance. Default is 1.
#' @param sigma_sq_0 Prior scale parameter for error variance. Default is 1.
#' @param n_MC Number of MCMC iterations (posterior samples). Default is 100.
#'
#' @return A list containing:
#'   \item{Lambda_samples}{Array of loading matrix posterior samples (genes × k × n_MC)}
#'   \item{sigma_sq_samples}{Vector of error variance posterior samples (length n_MC)}
#'   \item{rho}{Average correlation parameter used for MCMC tuning (scalar)}
#'
#' @export
#'
compute_posterior_samples_cc <- function(
    Y, Lambda_C, Lambda_N, tau_C, tau_N, sigma_sq, P_C,
    v0=1, sigma_sq_0=1, n_MC=100
){

  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(Lambda_C)

  Lambda <- Lambda_C + Lambda_N


  B <- compute_B(Lambda, sigma_sq)
  rho <- mean(B[lower.tri(B, diag = T)])
  print(paste('rho = ', rho))

  post_sample <- posterior_samples(
    Y, Lambda_C, Lambda_N, tau_C, tau_N, sigma_sq, P_C, v0, sigma_sq_0, n_MC,
    rho, rho)



  return(list(Lambda_samples = post_sample$Lambda_samples,
              sigma_sq_samples = post_sample$sigma_sq_samples,
              rho = rho))
}




predict_oos_Y <- function(Y_train, Y_test, impute_set, fit, P_C){

  # params posterior samples
  params_posterior_samples <- compute_posterior_samples_cc(
    Y_train, fit$Lambda_C, fit$Lambda_N, fit$tau_C,
    fit$tau_N, fit$sigma_sq, P_C
  )

  # impute latent factors
  M_posterior_samples <- sample_latent_factors(
    Y_test[,impute_set], params_posterior_samples$Lambda_samples[impute_set,,],
    params_posterior_samples$sigma_sq_samples
    )

  # get predictions
  predictions <- predict_Y_from_factors(
    M_posterior_samples$M_samples,
    params_posterior_samples$Lambda_samples[-impute_set,,],
    params_posterior_samples$sigma_sq_samples
  )
  return(list(Y_pred_mean = predictions$Y_mean,
              Y_pred_samples = predictions$Y_samples))
}

compute_covariance_posterior_mean <- function(Y, fit){
  n = dim(Y)[1]
  Lambda = fit$Lambda_C + fit$Lambda_N
  B <- compute_B(Lambda, fit$sigma_sq)
  rho <- mean(B[lower.tri(B, diag = T)])
  Lambda_outer_mean <- tcrossprod(Lambda) + fit$sigma_sq * rho *
    (1 / (n + 1/fit$tau_C) +1 / (n + 1/fit$tau_N)) * diag(nrow(Lambda))
  cov_mean <- Lambda_outer_mean + fit$sigma_sq * diag(nrow(Lambda))
  return(list(
    Lambda_outer_mean = Lambda_outer_mean, cov_mean = cov_mean
  ))
}


compute_Psiposterior = function(Lambda_samples, C){
  N_samples = dim(Lambda_samples)[3]
  p = dim(Lambda_samples)[1]
  k = dim(Lambda_samples)[2]
  tCCt = crossprod(C)
  s_cross_C = svd(tCCt)
  V_C = s_cross_C$u
  D_C = diag(as.vector(sqrt(s_cross_C$d)))
  D_C_inv = diag(as.vector(1/sqrt(s_cross_C$d)))
  U_C = C %*% V_C %*% D_C_inv
  P_C <- tcrossprod(U_C)
  Q_C <- diag(1, p, p) - P_C

  Psi_samples = array(0, dim=c(p, k, N_samples))
  for(i in 1:N_samples){
    Psi_samples[,,i] =  Q_C %*% Lambda_samples[,,i]
  }
  return(Psi_samples)
}


compute_Gammaposterior = function(Lambda_samples, C){
  N_samples = dim(Lambda_samples)[3]
  q = dim(C)[2]
  k = dim(Lambda_samples)[2]
  tCCt = crossprod(C)
  s_cross_C = svd(tCCt)
  V_C = s_cross_C$u
  D_C = diag(as.vector(sqrt(s_cross_C$d)))
  D_C_inv = diag(as.vector(1/sqrt(s_cross_C$d)))
  U_C = C %*% V_C %*% D_C_inv
  # (C'C)^{-1}C'
  C_left_pinv = V_C %*% D_C_inv %*% t(U_C)

  Gamma_samples = array(0, dim=c(q, k, N_samples))
  for(i in 1:N_samples){
    Gamma_samples[,,i] = C_left_pinv %*% Lambda_samples[,,i]
  }
  return(Gamma_samples)
}
