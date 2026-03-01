// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout << << std::endl; 
// For export to R: [[Rcpp::export]]

// ----------------------------------------------------------------------------
// Utilities.
// ----------------------------------------------------------------------------

// Check if value is in vector.
bool IsIn(double a, const arma::colvec& b) {
  return arma::any(b == a);
}

// Union of two time vectors; preserves order and uniqueness.
arma::colvec Union(const arma::colvec& a, arma::colvec b) {
  const double a0 = a(0);
  for (arma::uword i = 0; i < b.n_elem; ++i) {
    if (IsIn(b(i), a)) {
      b(i) = a0;
    }
  }
  return arma::unique(arma::join_cols(a, b));
}

// Truncate times at tau; optionally append tau if not present.
arma::colvec Truncate(const arma::colvec& time, double tau, bool add_tau = false) {
  arma::colvec unique_times = arma::unique(time);
  for (arma::uword i = 0; i < unique_times.n_elem; ++i) {
    if (unique_times(i) > tau) {
      unique_times(i) = tau;
    }
  }
  if (add_tau && !arma::any(unique_times == tau)) {
    const arma::uword n = unique_times.n_elem;
    arma::colvec out(n + 1);
    out.head(n) = unique_times;
    out(n) = tau;
    return arma::unique(out);
  }
  return arma::unique(unique_times);
}

// Add leading value; if already present, return x unchanged.
arma::colvec AddLeadVal(const arma::colvec& x, const double value) {
  if (arma::any(x == value)) {
    return x;
  }
  arma::colvec out(x.n_elem + 1);
  out(0) = value;
  out.tail(x.n_elem) = x;
  return out;
}


// ----------------------------------------------------------------------------
// Kaplan-Meier
// ----------------------------------------------------------------------------

// Tabulate KM: time, nar, surv, haz at each evaluation time.
arma::mat KaplanMeier(
    const arma::colvec& eval_times,
    const arma::colvec& status,
    const arma::colvec& time
) {
  const arma::uword n = time.n_elem;
  arma::colvec unique_times = Union(eval_times, arma::unique(time));
  const arma::uword n_unique_time = unique_times.n_elem;

  arma::colvec censor(n_unique_time);
  arma::colvec death(n_unique_time);
  arma::colvec nar(n_unique_time);
  double current_nar = static_cast<double>(n);

  for (arma::uword i = 0; i < n_unique_time; ++i) {
    const double current_time = unique_times(i);
    const arma::uvec idx = arma::find(time == current_time);
    nar(i) = current_nar;
    censor(i) = arma::accu(status.elem(idx) == 0.0);
    death(i) = arma::accu(status.elem(idx) == 1.0);
    current_nar -= censor(i) + death(i);
  }

  const arma::colvec haz = death / nar;
  const arma::colvec surv = arma::cumprod(1.0 - haz);

  const arma::uword n_eval = eval_times.n_elem;
  arma::colvec nar_out(n_eval);
  arma::colvec haz_out(n_eval);
  arma::colvec surv_out(n_eval);
  arma::uword ptr = 0;

  for (arma::uword i = 0; i < n_unique_time && ptr < n_eval; ++i) {
    if (IsIn(unique_times(i), eval_times)) {
      nar_out(ptr) = nar(i);
      haz_out(ptr) = haz(i);
      surv_out(ptr) = surv(i);
      ++ptr;
    }
  }

  return arma::join_rows(eval_times, nar_out, surv_out, haz_out);
}


// ----------------------------------------------------------------------------
// AUC
// ----------------------------------------------------------------------------

//' Calculate Restricted Mean Survival Time
//'
//' @param status Status, coded as 0 for censoring, 1 for death.
//' @param time Observation time.
//' @param extend Extend AUC calculation if tau exceeds max(time)?
//' @param tau Truncation time.
//' @return Numeric restricted mean survival time.
//' @export
// [[Rcpp::export]]
SEXP RMST(
  const arma::colvec& status,
  const arma::colvec& time,
  bool extend = false,
  Rcpp::Nullable<double> tau = R_NilValue
) {
  arma::colvec unique_times = arma::unique(time);
  double max_time = time.max();
  double trunc_time = max_time;
  if (tau.isNotNull()) {
    trunc_time = Rcpp::as<double>(tau);
    unique_times = Truncate(unique_times, trunc_time, false);
  }

  unique_times = AddLeadVal(unique_times, 0);
  const arma::uword n_times = unique_times.n_elem;

  const arma::mat km_mat = KaplanMeier(unique_times, status, time);
  const arma::colvec surv = km_mat.col(2);

  const arma::colvec delta_t = arma::diff(unique_times);
  const arma::colvec integrand = surv.head(n_times - 1);
  double auc = arma::sum(integrand % delta_t);

  if (extend && trunc_time > max_time) {
    auc += (trunc_time - max_time) * surv.min();
  }
  return Rcpp::wrap(auc);
}


// ----------------------------------------------------------------------------
// Martingales
// ----------------------------------------------------------------------------

// Subject (row) by evaluation time (col): dM[i,t] = dN[i,t] - Y[i,t]dA[t].
arma::mat CalcMartingale(
    const arma::colvec& eval_times,
    const arma::colvec& haz,
    const arma::colvec& status,
    const arma::colvec& time
) {
  const arma::uword n = time.n_elem;
  const arma::uword n_times = eval_times.n_elem;
  arma::mat dm = arma::zeros(n, n_times);

  for (arma::uword i = 0; i < n; ++i) {
    const double subj_time = time(i);
    const double subj_status = status(i);
    for (arma::uword j = 0; j < n_times; ++j) {
      if (eval_times(j) == subj_time && subj_status == 1.0) {
        dm(i, j) += 1.0;
      }
      if (eval_times(j) <= subj_time) {
        dm(i, j) -= haz(j);
      } else {
        break;
      }
    }
  }
  return dm;
}


// ----------------------------------------------------------------------------
// Integrate Kaplan Meier
// ----------------------------------------------------------------------------

// Area under the KM curve between a and b (step-function integral).
double IntegrateKM(
    double a,
    double b,
    const arma::colvec& eval_times,
    const arma::colvec& surv
) {
  if (a >= b) {
    return 0.0;
  }

  const arma::uvec idx_low = arma::find(eval_times <= a, 1, "last");
  const arma::uvec idx_high = arma::find(eval_times <= b, 1, "last");
  const double lower_surv = idx_low.n_elem > 0 ? surv(idx_low(0)) : 1.0;
  const double upper_surv = idx_high.n_elem > 0 ? surv(idx_high(0)) : 0.0;

  const arma::uvec mid = arma::find((eval_times > a) % (eval_times < b));
  arma::colvec int_times = arma::join_cols(
    eval_times.elem(mid),
    arma::vec({a, b})
  );
  int_times = arma::unique(int_times);
  const arma::uword n_t = int_times.n_elem;
  arma::colvec int_vals(n_t);

  for (arma::uword i = 0; i < n_t; ++i) {
    const arma::uvec k = arma::find(eval_times == int_times(i));
    if (i == 0) {
      int_vals(i) = lower_surv;
    } else if (i == n_t - 1) {
      int_vals(i) = upper_surv;
    } else {
      int_vals(i) = k.n_elem > 0 ? surv(k(0)) : lower_surv;
    }
  }

  const arma::colvec dt = arma::diff(int_times);
  return arma::as_scalar(arma::sum(int_vals.head(dt.n_elem) % dt));
}


// ----------------------------------------------------------------------------
// RMST influence function
// ----------------------------------------------------------------------------

//' RMST Influence Function
//'
//' Influence function of the restricted mean survival time at time t.
//'
//' @param status Status, coded as 0 for censoring, 1 for death.
//' @param time Observation time.
//' @param trunc_time Truncation time.
//' @return Numeric vector of influence function values for each observation.
// [[Rcpp::export]]
SEXP CalcPsiRMST(
  const arma::colvec& status,
  const arma::colvec& time,
  double trunc_time
) {
  arma::colvec unique_times = AddLeadVal(arma::unique(time), 0.0);
  unique_times = Truncate(unique_times, trunc_time);
  const arma::uword n_subj = time.n_elem;

  const arma::mat km_mat = KaplanMeier(unique_times, status, time);
  const arma::colvec km_times = km_mat.col(0);
  const arma::colvec nar = km_mat.col(1);
  const arma::colvec surv = km_mat.col(2);
  const arma::colvec haz = km_mat.col(3);
  const arma::colvec par = nar / static_cast<double>(n_subj);

  arma::mat mart = CalcMartingale(km_times, haz, status, time).t();
  const arma::uword n_times = km_times.n_elem;
  arma::colvec mu(n_times);
  for (arma::uword j = 0; j < n_times; ++j) {
    mu(j) = IntegrateKM(km_times(j), trunc_time, km_times, surv);
  }

  arma::colvec psi(n_subj);
  for (arma::uword i = 0; i < n_subj; ++i) {
    double sum_val = 0.0;
    for (arma::uword j = 0; j < n_times; ++j) {
      if (par(j) > 0.0) {
        sum_val += (mu(j) / par(j)) * mart(j, i);
      }
    }
    psi(i) = -sum_val;
  }
  return Rcpp::wrap(psi);
}


// ----------------------------------------------------------------------------
// Run perturbations.
// ----------------------------------------------------------------------------

//' Generate Perturbations
//'
//' @param psi Per-subject influence values.
//' @param n_boot Number of perturbations.
//' @return Numeric vector.
//' @section Notes:
//' The random seed should be set in R prior to calling this function.
// [[Rcpp::export]]
SEXP GenPerturb(
  const arma::colvec& psi,
  const int n_boot
) {
  const arma::uword n = psi.n_elem;
  arma::mat W = arma::randn(n, n_boot);
  arma::colvec out(n_boot);
  for (int b = 0; b < n_boot; ++b) {
    out(b) = arma::mean(psi % W.col(b));
  }
  return Rcpp::wrap(out);
}
