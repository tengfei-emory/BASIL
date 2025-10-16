
#' Estimate the latent dimension using JIC
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



#' Compute point estimates for BASIL model parameters
#'
#' @param Y A numeric matrix of gene expression data (samples × genes)
#' @param C A binary gene membership matrix (genes × pathways).
#' @param k Integer specifying the number of latent factors. If NA (default),
#'   the function will estimate k using \code{estimate_latent_dimension}.
#' @param k_max Maximum number of factors to consider when estimating k.
#'   Only used if k = NA. Default is 50.
#' @param v0 Prior hyperparameter for the error variance. Default is 1.
#' @param sigma_sq0 Prior hyperparameter for the error variance. Default is 1.
#'
#' @return A list containing:
#'   \item{M}{Matrix of latent factor scores (samples × k)}
#'   \item{Lambda_C}{Loading matrix in the columns space of C, i.e. pathway-guided loadings (genes × k)}
#'   \item{Lambda_N}{Loading matrix in the null space of C, i.e. unstructured loadings (genes × k)}
#'   \item{sigma_sq}{Estimated error variance (scalar)}
#'   \item{tau_gamma}{Prior variance parameter for Lambda_C (scalar)}
#'   \item{tau_psi}{Prior  variance parameter for Lambda_N (scalar)}
#'   \item{k}{Number of latent factors used (integer)}
#'   \item{P_C}{Projection matrix onto the column space of C (genes × genes)}
#'
#' @export
#'
BASIL_point_estimates <- function(Y, C, k = NA, k_max = 50, v0 = 1, sigma_sq0 = 1){

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
  tau_gamma <- sum((P_V_D)^2 ) / (k * sigma_sq_hat * q) / n
  print(paste('tau_gamma  = ', tau_gamma))
  tau_psi <- sum((Q_V_D)^2 ) / (k * sigma_sq_hat * (p-q))/ n
  print(paste('tau_psi  = ', tau_psi))

  # point estimates
  M <- sqrt(n) * U
  Lambda_C <- (P_V_D) * sqrt(n) / (n + 1/tau_gamma)
  Lambda_N <- (Q_V_D) * sqrt(n) / (n + 1/tau_psi)
  sigma_sq <- (v0 * sigma_sq0 + sum((Y-tcrossprod(M, Lambda_C + Lambda_N) )^2) ) /
    (v0 + n*p -2)

  return(list(M=M, Lambda_C=Lambda_C, Lambda_N=Lambda_N, sigma_sq=sigma_sq,
              tau_gamma=tau_gamma, tau_psi=tau_psi, k=k, P_C=P_C))
}



#' Compute posterior samples for BASIL model parameters with coverage correction
#'
#' @param Y A numeric matrix of gene expression data (samples × genes)
#' @param Lambda_C Loading matrix in the columns space of C, i.e. pathway-guided loadings (genes × k)
#' @param Lambda_N Loading matrix in the null space of C, i.e. unstructured loadings (genes × k)
#' @param tau_gamma Prior variance parameter for Lambda_C (scalar)
#' @param tau_psi Prior variance parameter for Lambda_N (scalar)
#' @param sigma_sq Error variance estimate
#' @param P_C Projection matrix onto the column space of C (genes × genes)
#' @param v0 Prior hyperparameter for the error variance. Default is 1.
#' @param sigma_sq_0 Prior hyperparameter for the error variance. Default is 1.
#' @param n_MC Number of MCMC iterations (posterior samples). Default is 100.
#'
#' @return A list containing:
#'   \item{Lambda_samples}{Array of loading matrix posterior samples (genes × k × n_MC)}
#'   \item{sigma_sq_samples}{Vector of error variance posterior samples (length n_MC)}
#'   \item{rho}{Coverage correction term (scalar)}
#'
#' @export
#'
compute_posterior_samples_cc <- function(Y, Lambda_C, Lambda_N, tau_gamma, tau_psi,
                                         sigma_sq, P_C, v0=1, sigma_sq_0=1, n_MC=100){

  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(Lambda_C)

  Lambda <- Lambda_C + Lambda_N
  B <- compute_B(Lambda, sigma_sq)
  rho <- mean(B[lower.tri(B, diag = T)])
  print(paste('rho = ', rho))

  post_sample <- posterior_samples(
    Y, Lambda_C, Lambda_N, tau_gamma, tau_psi, sigma_sq, P_C, v0, sigma_sq_0, n_MC,
    rho, rho)

  return(list(Lambda_samples = post_sample$Lambda_samples,
              sigma_sq_samples = post_sample$sigma_sq_samples,
              rho = rho))
}



#' Compute posterior mean of the covariance matrix
#'
#' @param Y A numeric matrix of gene expression data (samples × genes)
#' @param fit A list containing model fit from \code{compute_point_estimates},
#'   with components:
#'   \describe{
#'     \item{Lambda_C}{Pathway-guided loading matrix (genes × k)}
#'     \item{Lambda_N}{Non-pathway loading matrix (genes × k)}
#'     \item{sigma_sq}{Error variance estimate}
#'     \item{tau_gamma}{Prior precision for pathway-guided loadings}
#'     \item{tau_psi}{Prior precision for non-pathway loadings}
#'   }
#'
#' @return A list containing:
#'   \item{Lambda_outer_mean}{Posterior mean of the low-rank covariance
#'     component Lambda * Lambda^T}
#'   \item{cov_mean}{Posterior mean of the full gene-gene covariance matrix}
#'
#' @export
#'
compute_covariance_posterior_mean <- function(Y, fit){
  n <- dim(Y)[1]
  Lambda <- fit$Lambda_C + fit$Lambda_N
  B <- compute_B(Lambda, fit$sigma_sq)
  rho <- mean(B[lower.tri(B, diag = T)])
  Lambda_outer_mean <- tcrossprod(Lambda) + fit$sigma_sq * rho *
    (1 / (n + 1/fit$tau_gamma) +1 / (n + 1/fit$tau_psi)) * diag(nrow(Lambda))
  cov_mean <- Lambda_outer_mean + fit$sigma_sq * diag(nrow(Lambda))

  return(list(Lambda_outer_mean = Lambda_outer_mean, cov_mean = cov_mean))
}


