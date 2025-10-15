#include <RcppArmadillo.h>
#include <Rcpp.h>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//' Compute posterior samples of the gene gene correlation matrix
//'
//' @param Lambda_samples Array of loading matrix posterior samples (genes × k × n_samples).
//' @param sigma_sq_samples Vector of error variance posterior samples (length n_samples).
//' @param samples Logical flag indicating whether to return all posterior samples.
//'
//' @return A list containing:
//'   \item{posterior_mean}{Posterior mean of the correlation matrix (genes × genes).}
//'   \item{posterior_samples}{(Only if samples = TRUE) Array of correlation matrix
//'     posterior samples (genes × genes × n_samples).}
//'
//' @export
//'
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


//' Compute posterior samples of Gamma
//'
//' @param Lambda_samples Array of posterior samples of the loading matrix
//'   (genes × k × n_samples)
//' @param C A binary gene membership matrix (genes × pathways)
//' @param R_in Optional rotation matrix (k × k) to apply to Lambda before
//'   computing Gamma. If NULL, no rotation is applied. Default is NULL.
//'
//' @return An array of posterior samples of Gamma (pathways × k × n_samples).
//'
//' @export
//'
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


//' Compute posterior samples of Psi
//'
//' @param Lambda_samples Array of posterior samples of the loading matrix
//'   (genes × k × n_samples)
//' @param P_C Projection matrix onto the column space of C (genes × genes).
//'    Typically obtained from \code{compute_point_estimates}.
//' @param R_in Optional rotation matrix (k × k) to apply to Lambda before
//'   computing Psi. If not provided or NULL (default), no rotation is applied.
//'
//' @return An array of posterior samples of Psi (genes × k × n_samples),
//'
//' @export
//'
// [[Rcpp::export]]
arma::cube compute_Psi_samples_cpp(
   const arma::cube& Lambda_samples,
   const arma::mat& P_C,
   Rcpp::Nullable<arma::mat> R_in = R_NilValue
) {
  const arma::uword p    = Lambda_samples.n_rows;
  const arma::uword k    = Lambda_samples.n_cols;
  const arma::uword n_MC = Lambda_samples.n_slices;

  // Check dimensions
  if (P_C.n_rows != p || P_C.n_cols != p) {
    Rcpp::stop("P_C must be a square matrix with dimensions matching nrow(Lambda_samples).");
  }

  // Compute Q_C = I - P_C (projection onto null space)
  arma::mat Q_C = arma::eye<arma::mat>(p, p) - P_C;

  // Check if rotation matrix is provided
  arma::mat Rmat;
  bool use_R = false;
  if (R_in.isNotNull()) {
    Rmat = Rcpp::as<arma::mat>(R_in);
    if (Rmat.n_rows != k || Rmat.n_cols != k) {
      Rcpp::stop("R_in must be a square matrix with dimensions matching ncol(Lambda_samples).");
    }
    use_R = true;
  }

  // Initialize output
  arma::cube Psi_samples(p, k, n_MC, arma::fill::none);

  // Compute Psi for each sample
  for (arma::uword t = 0; t < n_MC; ++t) {
    arma::mat L = Lambda_samples.slice(t);
    if (use_R) {
      L = L * Rmat;
    }
    Psi_samples.slice(t) = Q_C * L;
  }

  return Psi_samples;
}


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







