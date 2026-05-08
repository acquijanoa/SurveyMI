#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

namespace {
arma::mat safe_inv(const arma::mat& m) {
  arma::mat out;
  bool ok = arma::inv_sympd(out, m);
  if (!ok) {
    out = arma::pinv(m);
  }
  return out;
}

arma::mat ridge_eye(std::size_t n, double raw_scale) {
  const double jitter =
      raw_scale <= 0.0
          ? 1e-10
          : (1e-10 * raw_scale); // stabilises XTWX near rank-deficient crosses
  return jitter * arma::eye<arma::mat>(static_cast<unsigned int>(n),
                                        static_cast<unsigned int>(n));
}

arma::vec sys_solve_stable(const arma::mat& xtwx_raw, const arma::vec& rhs) {
  const double dn = static_cast<double>(std::max(1u, xtwx_raw.n_rows));
  const double tr = arma::norm(arma::diagvec(xtwx_raw), 1) / dn;
  arma::mat a = xtwx_raw + ridge_eye(xtwx_raw.n_rows, tr);
  arma::vec coef;
  const bool ok = arma::solve(
      coef,
      a,
      rhs,
      arma::solve_opts::fast + arma::solve_opts::likely_sympd
  );
  if (!ok || !coef.is_finite()) {
    coef = arma::pinv(a) * rhs;
  }
  return coef;
}
} // namespace

// [[Rcpp::export]]
Rcpp::List weighted_gaussian_cpp(const arma::mat& x,
                                 const arma::vec& y,
                                 const arma::vec& w) {
  const arma::vec sw = arma::sqrt(arma::clamp(w, 0.0, arma::datum::inf));
  const arma::mat xw = x.each_col() % sw;
  const arma::vec yw = y % sw;

  const arma::mat xtwx = xw.t() * xw;
  const arma::vec xtwy = xw.t() * yw;
  const arma::vec beta = sys_solve_stable(xtwx, xtwy);

  const arma::vec resid = y - x * beta;
  const double rss = arma::dot(w, resid % resid);
  const double n_eff = arma::accu(w > 0.0);
  const double p = x.n_cols;
  const double sigma2 = (n_eff > p) ? (rss / (n_eff - p)) : std::max(rss, 1e-8);

  const arma::mat vcov = sigma2 * safe_inv(xtwx);
  return Rcpp::List::create(
    Rcpp::Named("coef") = beta,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("vcov") = vcov
  );
}

// [[Rcpp::export]]
Rcpp::List weighted_binomial_cpp(const arma::mat& x,
                                 const arma::vec& y,
                                 const arma::vec& w,
                                 const int maxit = 50,
                                 const double tol = 1e-8) {
  arma::vec beta(x.n_cols, arma::fill::zeros);
  arma::mat xtwx;
  arma::vec mu;

  for (int it = 0; it < maxit; ++it) {
    const arma::vec eta = x * beta;
    mu = 1.0 / (1.0 + arma::exp(-eta));
    mu = arma::clamp(mu, 1e-8, 1.0 - 1e-8);

    const arma::vec work_w = w % (mu % (1.0 - mu));
    const arma::vec z = eta + (y - mu) / (mu % (1.0 - mu));

    const arma::vec sw = arma::sqrt(arma::clamp(work_w, 0.0, arma::datum::inf));
    const arma::mat xw = x.each_col() % sw;
    const arma::vec zw = z % sw;

    xtwx = xw.t() * xw;
    const arma::vec xtwz = xw.t() * zw;
    const arma::vec beta_new = sys_solve_stable(xtwx, xtwz);

    if (arma::norm(beta_new - beta, 2) < tol) {
      beta = beta_new;
      break;
    }
    beta = beta_new;
  }

  const arma::mat vcov = safe_inv(xtwx);
  return Rcpp::List::create(
    Rcpp::Named("coef") = beta,
    Rcpp::Named("vcov") = vcov
  );
}

// [[Rcpp::export]]
Rcpp::List weighted_poisson_cpp(const arma::mat& x,
                                const arma::vec& y,
                                const arma::vec& w,
                                const int maxit = 50,
                                const double tol = 1e-8) {
  arma::vec beta(x.n_cols, arma::fill::zeros);
  arma::mat xtwx;
  arma::vec mu;

  for (int it = 0; it < maxit; ++it) {
    const arma::vec eta = x * beta;
    mu = arma::exp(eta);
    mu = arma::clamp(mu, 1e-8, arma::datum::inf);

    const arma::vec work_w = w % mu;
    const arma::vec z = eta + (y - mu) / mu;

    const arma::vec sw = arma::sqrt(arma::clamp(work_w, 0.0, arma::datum::inf));
    const arma::mat xw = x.each_col() % sw;
    const arma::vec zw = z % sw;

    xtwx = xw.t() * xw;
    const arma::vec xtwz = xw.t() * zw;
    const arma::vec beta_new = sys_solve_stable(xtwx, xtwz);

    if (arma::norm(beta_new - beta, 2) < tol) {
      beta = beta_new;
      break;
    }
    beta = beta_new;
  }

  const arma::mat vcov = safe_inv(xtwx);
  return Rcpp::List::create(
    Rcpp::Named("coef") = beta,
    Rcpp::Named("vcov") = vcov
  );
}

namespace {

// 0 = observed, 1 = missing (NaN)
std::string get_missing_pattern(const arma::rowvec& x) {
  std::string pattern;
  pattern.reserve(x.n_elem);
  for (unsigned int j = 0; j < x.n_elem; j++) {
    pattern += std::isnan(x(j)) ? '1' : '0';
  }
  return pattern;
}

} // namespace

// One weighted EM step for MVN with arbitrary missingness patterns
// [[Rcpp::export]]
Rcpp::List run_fast_mvn_em_step(arma::mat data,
                                const arma::vec& mu,
                                const arma::mat& Sigma,
                                const arma::vec& weights) {
  const int n = static_cast<int>(data.n_rows);
  const int p = static_cast<int>(data.n_cols);

  arma::mat data_imp = data;
  arma::mat expected_cov = arma::zeros(p, p);
  arma::vec expected_mean = arma::zeros(p);
  double sum_w = 0.0;

  std::unordered_map<std::string, std::vector<int>> pattern_groups;
  if (n > 0) {
    pattern_groups.reserve(static_cast<unsigned long>(n));
  }

  for (int i = 0; i < n; i++) {
    if (weights(static_cast<unsigned int>(i)) <= 0.0) {
      continue;
    }
    std::string pat =
        get_missing_pattern(data.row(static_cast<unsigned int>(i)));
    pattern_groups[std::move(pat)].push_back(i);
  }

  for (const auto& pair : pattern_groups) {
    const std::string& pat = pair.first;
    const std::vector<int>& row_indices = pair.second;

    std::vector<arma::uword> obs_vec;
    std::vector<arma::uword> mis_vec;
    obs_vec.reserve(static_cast<unsigned long>(pat.size()));
    mis_vec.reserve(static_cast<unsigned long>(pat.size()));
    for (unsigned int j = 0; j < pat.size(); j++) {
      if (pat[j] == '0') {
        obs_vec.push_back(j);
      } else {
        mis_vec.push_back(j);
      }
    }
    arma::uvec obs = arma::conv_to<arma::uvec>::from(obs_vec);
    arma::uvec mis = arma::conv_to<arma::uvec>::from(mis_vec);

    if (mis.n_elem == 0) {
      for (int i : row_indices) {
        const unsigned int ui = static_cast<unsigned int>(i);
        arma::vec x_i = data.row(ui).t();
        const double wi = weights(ui);
        expected_mean += wi * x_i;
        expected_cov += wi * (x_i * x_i.t());
        sum_w += wi;
      }
      continue;
    }

    if (obs.n_elem == 0) {
      for (int i : row_indices) {
        const unsigned int ui = static_cast<unsigned int>(i);
        data_imp.row(ui) = mu.t();
        const double wi = weights(ui);
        expected_mean += wi * mu;
        expected_cov += wi * (mu * mu.t() + Sigma);
        sum_w += wi;
      }
      continue;
    }

    const arma::vec mu_obs = mu.elem(obs);
    const arma::vec mu_mis = mu.elem(mis);

    const arma::mat S_obs_obs = Sigma.submat(obs, obs);
    const arma::mat S_mis_obs = Sigma.submat(mis, obs);
    const arma::mat S_mis_mis = Sigma.submat(mis, mis);

    arma::mat S_obs_obs_inv;
    const bool inv_ok = arma::inv_sympd(S_obs_obs_inv, S_obs_obs);
    if (!inv_ok) {
      S_obs_obs_inv = arma::pinv(S_obs_obs);
    }

    const arma::mat projection = S_mis_obs * S_obs_obs_inv;
    const arma::mat V_mis = S_mis_mis - projection * S_mis_obs.t();

    arma::mat V_i = arma::zeros(p, p);
    V_i.submat(mis, mis) = V_mis;

    for (int i : row_indices) {
      const unsigned int ui = static_cast<unsigned int>(i);
      arma::vec x_i = data.row(ui).t();
      const arma::vec x_obs = x_i.elem(obs);
      const arma::vec x_mis_hat = mu_mis + projection * (x_obs - mu_obs);
      x_i.elem(mis) = x_mis_hat;
      data_imp.row(ui) = x_i.t();
      const double wi = weights(ui);
      expected_mean += wi * x_i;
      expected_cov += wi * (x_i * x_i.t() + V_i);
      sum_w += wi;
    }
  }

  if (sum_w <= 0.0) {
    return Rcpp::List::create(Rcpp::Named("mu") = mu,
                              Rcpp::Named("Sigma") = Sigma,
                              Rcpp::Named("data_imp") = data_imp,
                              Rcpp::Named("sum_w") = sum_w);
  }

  const arma::vec new_mu = expected_mean / sum_w;
  arma::mat new_Sigma = (expected_cov / sum_w) - (new_mu * new_mu.t());
  new_Sigma.diag() += 1e-6;

  return Rcpp::List::create(Rcpp::Named("mu") = new_mu,
                            Rcpp::Named("Sigma") = new_Sigma,
                            Rcpp::Named("data_imp") = data_imp,
                            Rcpp::Named("sum_w") = sum_w);
}
