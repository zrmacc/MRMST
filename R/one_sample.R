# Purpose: One-sample multiple RMST calculation.
# Updated: 2023-09-29


#' One Sample Calculation
#'
#' @param statuses Array with 1 or more status variables.
#' @param times Array with 1 or more event-time variables.
#' @param alpha Type I error.
#' @param extend Extend the Kaplan-Meier curves if `tau` is set beyond the
#'   maximum observation time?
#' @param n_boot Number of perturbations for variance calculation.
#' @param seed Optional integer. If set, the RNG state is fixed before running
#'   perturbations so that results are reproducible.
#' @param tau Truncation time.
#' @param weights Optional weight vector. Defaults to equal weights.
#' @export
OneSample <- function(
  statuses,
  times,
  alpha = 0.05,
  extend = FALSE,
  n_boot = 2000,
  seed = NULL,
  tau = NULL,
  weights = NULL
) {
  
  statuses <- as.matrix(statuses)
  times <- as.matrix(times)
  n_comp <- ncol(times)

  InputsChecks(statuses = statuses, times = times)

  if (is.null(weights)) {
    weights <- rep(1, n_comp)
  }
  if (is.null(tau)) {
    tau <- MaxEventTime(statuses = statuses, times = times)
  }
  InputsChecks(
    statuses = statuses,
    times = times,
    tau = tau,
    weights = weights,
    n_comp = n_comp
  )

  n_subj <- nrow(times)
  aucs <- numeric(n_comp)
  psi_mat <- matrix(0, nrow = n_subj, ncol = n_comp)
  for (i in seq_len(n_comp)) {
    aucs[i] <- RMST(
      status = statuses[, i],
      time = times[, i],
      extend = extend,
      tau = tau
    )
    psi_mat[, i] <- CalcPsiRMST(
      status = statuses[, i],
      time = times[, i],
      trunc_time = tau
    )
  }
  obs_stat <- sum(weights * aucs)
  psi <- as.numeric(psi_mat %*% weights)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  perturb <- GenPerturb(psi = psi, n_boot = n_boot)
  
  # Output.
  z <- stats::qnorm(1 - alpha / 2)
  table <- data.frame(
    n = nrow(times),
    tau = tau,
    k = n_comp,
    auc = obs_stat,
    se = stats::sd(perturb)
  )
  table$lower <- table$auc - z * table$se
  table$upper <- table$auc + z * table$se
  
  out <- methods::new(
    "OneSample",
    AUC = table$auc,
    K = n_comp,
    SE = table$se,
    Table = table,
    Tau = tau
  )
  return(out)
}

