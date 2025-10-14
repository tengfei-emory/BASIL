#include <RcppArmadillo.h>
#include <Rcpp.h>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat compute_B(const arma::mat& Lambda, const double& sigma_sq) {
  arma::mat Lambda_outer = Lambda * Lambda.t(); 
  int p = Lambda.n_rows; 
  arma::mat B = arma::zeros<arma::mat>(p, p);
  
  for (int j = 0; j < p - 1; ++j) {
    for (int l = j + 1; l < p; ++l) {
      B(j, l) = sqrt(1 + (Lambda_outer(j, j) * Lambda_outer(l, l) + Lambda_outer(j, l) * Lambda_outer(j, l)) / 
        (sigma_sq * Lambda_outer(l, l) + sigma_sq * Lambda_outer(j, j)));
      B(l, j) = B(j, l); 
    }
    B(j, j) = sqrt(1 + Lambda_outer(j, j) / (2 * sigma_sq));
  }
  B(p - 1, p - 1) = sqrt(1 + Lambda_outer(p - 1, p - 1) / (2 * sigma_sq));
  
  return B;
}

// [[Rcpp::export]]
List posterior_samples(const arma::mat& Y, 
                   const arma::mat& Lambda_C_hat,
                   const arma::mat& Lambda_N_hat,
                   double tau_C,
                   double tau_N,
                   double sigma_sq_hat, 
                   const arma::mat P_C,
                   double v_0 = 1, 
                   double sigma_sq_0 = 1, 
                   int N_mc = 100, 
                   double rho_C = 1, 
                   double rho_N = 1) {
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  int k = Lambda_C_hat.n_cols;

  arma::cube Lambda_samples(p, k, N_mc, arma::fill::zeros);

  arma::mat Lambda_t(p, k, arma::fill::zeros);
  arma::vec sigma_sq_samples(N_mc, arma::fill::zeros);
  double residuals;
  
  arma::mat Q_C = arma::eye(p, p) - P_C;
  
  double SST = arma::accu(arma::square(Y));
  double SSE = (n + 1/tau_C) * arma::accu(arma::square(Lambda_C_hat)) + 
    (n + 1/tau_N) * arma::accu(arma::square(Lambda_N_hat));
  double v_n = v_0 + n*p;
  double gamma_n = v_0 * sigma_sq_0 + SST - SSE;

  //arma::mat D = arma::diagmat((residuals + gamma_0 * delta_0) / (gamma_0 + n - 2)) * k * rho * rho / (n + 1 / tau_2);
  //Lambda_outer_mean = mu_js * mu_js.t() + D;
  
  // posterior sampling
  for (int it = 0; it < N_mc; ++it) {
    sigma_sq_samples(it) = 1.0 / R::rgamma(v_n / 2, 2.0 / gamma_n);
    
    arma::mat eps_C = sqrt(sigma_sq_samples(it)) * rho_C * arma::randn<arma::mat>(p, k) / sqrt(n + 1/tau_C);
    arma::mat eps_C_P = P_C * eps_C;
    arma::mat eps_N = sqrt(sigma_sq_samples(it)) * rho_N * arma::randn<arma::mat>(p, k) / sqrt(n + 1/tau_N);
    arma::mat eps_N_Q = Q_C * eps_N;
    
    Lambda_t = Lambda_C_hat + Lambda_N_hat + eps_C_P + eps_N_Q;
    
    Lambda_samples.slice(it) = Lambda_t;
  }
  
  return List::create(
    Named("Lambda_samples") = Lambda_samples,
    Named("sigma_sq_samples") = sigma_sq_samples
  );
}


// [[Rcpp::export]]
arma::cube sample_Lambda_outer(arma::cube Lambda_samples){
  
  int n_MC = Lambda_samples.n_slices;
  int k = Lambda_samples.n_cols;
  int p = Lambda_samples.n_rows;
  
  arma::mat Lambda(p, k);
  arma::mat Lambda_outer(p, p);
  arma::cube Lambda_outer_samples(p, p, n_MC);
  
  for(int s=0; s<n_MC; ++s) {
    Lambda = Lambda_samples.slice(s);
    Lambda_outer = Lambda * Lambda.t();
    Lambda_outer_samples.slice(s) = Lambda_outer;
  }
  
  return Lambda_outer_samples;
}



// [[Rcpp::export]]
List sample_latent_factors(
    arma::mat Y, arma::cube Lambda_samples, arma::vec sigma_sq_samples){
  
  int n_MC = Lambda_samples.n_slices;
  int k = Lambda_samples.n_cols;
  int p = Lambda_samples.n_rows;
  int n = Y.n_rows;
  
  arma::mat M_mean(n, k);
  arma::cube M_samples(n, k, n_MC);
  
  double sigma_sq;
  arma::mat Lambda(p, k);
  arma::mat prec_mat(k, k);
  
  arma::mat U;
  arma::vec d;
  arma::mat V;
  
  arma::vec inv_d;
  arma::vec sqrt_inv_d;
  
  arma::mat var_mat(k,k);
  arma::mat var_sqrt(k,k);
  
  arma::mat mean_t(n, k);
  arma::mat noise(n, k);
  
  for(int s=0; s<n_MC; ++s) {
    
    sigma_sq = sigma_sq_samples(s);
    Lambda = Lambda_samples.slice(s);
    prec_mat = arma::eye(k,k) + 1/sigma_sq * (Lambda.t() * Lambda);
    
    svd(U, d, V, prec_mat);
    
    // var_mat = U * diag(1/d) * U^T
    inv_d = 1.0 / d;
    var_mat = U * diagmat(inv_d) * U.t();
    
    // var_sqrt = U * diag(sqrt(1/d)) * U^T
    sqrt_inv_d = sqrt(1.0 / d);
    var_sqrt = U * diagmat(sqrt_inv_d) * U.t();
    
    // mean_t = (1 / sigma_sq) * Y * Lambda * var_mat
    mean_t = (1.0 / sigma_sq) * Y * Lambda * var_mat;
    
    // Update M_mean
    M_mean += mean_t;
    
    // Generate Gaussian noise
    noise = arma::randn(n, k);
    
    // Add noise scaled by var_sq and store in cube
    M_samples.slice(s) = mean_t + noise * var_sqrt;
    
  }
  
  M_mean = M_mean / n_MC;
  
  return List::create(
    Named("M_mean") = M_mean,
    Named("M_samples") = M_samples
  );
  
}

// [[Rcpp::export]]
List predict_Y_from_factors(arma::cube M_samples, arma::cube Lambda_samples, arma::vec sigma_sq_samples){
  
  int n_MC = Lambda_samples.n_slices;
  int k = Lambda_samples.n_cols;
  int p = Lambda_samples.n_rows;
  int n = M_samples.n_rows;
  
  arma::mat Y_mean(n, p);
  arma::mat Y_mean_s(n, p);
  arma::mat Y_sample_s(n, p);
  arma::cube Y_samples(n, p, n_MC);
  
  for(int s=0; s<n_MC; ++s) {
    Y_mean_s = M_samples.slice(s) * (Lambda_samples.slice(s)).t();
    Y_sample_s = Y_mean_s + sqrt(sigma_sq_samples(s)) * arma::randn(n, p);
    Y_samples.slice(s) = Y_sample_s;
    Y_mean = Y_mean + Y_mean_s;
  }
  Y_mean = Y_mean / n_MC;
  
  return List::create(
    Named("Y_mean") = Y_mean,
    Named("Y_samples") = Y_samples
  );
}


// [[Rcpp::export]]
List compute_covariance_posterior_samples_cc(
    const arma::cube& Lambda_samples, const arma::vec& sigma_sq_samples,
    bool samples = true, bool only_low_rank = false) {
  int p     = Lambda_samples.n_rows;
  int n_MC  = Lambda_samples.n_slices;
  
  arma::mat cov_mean(p, p, fill::zeros);
  arma::cube cov_samples;                         
  if (samples) cov_samples.set_size(p, p, n_MC);
  
  for (int t = 0; t < n_MC; ++t) {
    const arma::mat& L_t = Lambda_samples.slice(t);     
    arma::mat cov_t = L_t * L_t.t();                   
    if (!only_low_rank) {
      cov_t.diag() += sigma_sq_samples[t];
    }    
    cov_mean += cov_t / static_cast<double>(n_MC);
    if (samples)  {
      cov_samples.slice(t) = std::move(cov_t);
    }
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["posterior_mean"] = cov_mean
  );
  if (samples) {
    out["posterior_samples"] = cov_samples;             
  }
  return out;
}


// [[Rcpp::export]]
List compute_correlation_posterior_samples_cc(
    const arma::cube& Lambda_samples,
    const arma::vec&  sigma_sq_samples,
    bool samples = true
) {
  const int p    = Lambda_samples.n_rows;
  const int n_MC = Lambda_samples.n_slices;
  
  arma::mat  corr_mean(p, p, arma::fill::zeros);
  arma::cube corr_samples;
  if (samples) {
    corr_samples.set_size(p, p, n_MC);
  }
  
  const double eps = 1e-12;
  
  auto cov_to_corr = [&](const arma::mat& S) -> arma::mat {
    arma::vec sd = arma::sqrt(S.diag());
    sd.transform([&](double v){ return (v > 0.0) ? v : eps; });
    arma::mat C = S;
    C.each_col() /= sd;
    C.each_row() /= sd.t();
    C.diag().ones();
    return C;
  };
  
  for (int t = 0; t < n_MC; ++t) {
    const arma::mat& L_t = Lambda_samples.slice(t);
    arma::mat cov_t = L_t * L_t.t();
    cov_t.diag() += sigma_sq_samples[t];
    arma::mat corr_t = cov_to_corr(cov_t);
    
    corr_mean += corr_t / static_cast<double>(n_MC);
    if (samples) {
      corr_samples.slice(t) = std::move(corr_t);
    }
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["posterior_mean"] = corr_mean
  );
  if (samples) out["posterior_samples"] = corr_samples;
  return out;
}

// [[Rcpp::export]]
arma::cube compute_Gamma_samples_cpp(
    const arma::cube& Lambda_samples, const arma::mat& C,
    Rcpp::Nullable<arma::mat> R_in = R_NilValue
) {
  
  arma::mat U, V;
  arma::vec s;
  if (!arma::svd(U, s, V, C)) {
    stop("SVD failed on C.");
  }
  
  const double tol = 1e-5;
  arma::uword d = arma::sum(s > tol);

  arma::mat Usub = U.cols(0, d - 1);
  arma::mat Vsub = V.cols(0, d - 1);
  arma::vec sinv = 1.0 / s.head(d);
  arma::mat Dinv = arma::diagmat(sinv);
  arma::mat ols = Vsub * Dinv * Usub.t();
  
  const arma::uword p    = Lambda_samples.n_rows;
  const arma::uword k    = Lambda_samples.n_cols;
  const arma::uword n_MC = Lambda_samples.n_slices;
  const arma::uword q    = C.n_cols;
  
  if (C.n_rows != p) stop("nrow(C) must equal nrow(Lambda_samples).");
  
  arma::mat Rmat;
  bool use_R = false;
  if (R_in.isNotNull()) {
    Rmat = as<arma::mat>(R_in);
    if (Rmat.n_rows != k) stop("nrow(R) must equal ncol(Lambda_samples).");
    if (Rmat.n_cols != k) stop("ncol(R) must equal ncol(Lambda_samples).");
    use_R = true;
  }
  
  arma::cube Gamma_samples(q, k, n_MC, arma::fill::none);
  
  for (arma::uword t = 0; t < n_MC; ++t) {
    arma::mat L = Lambda_samples.slice(t);
    if (use_R) L = L * Rmat;
    Gamma_samples.slice(t) = ols * L;
  }
  
  return Gamma_samples;
}



