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
#' @param tau Truncation time.
#' @param weights Optional weight vector. Defaults to equal weights.
#' @export
OneSample <- function(
  statuses,
  times,
  alpha = 0.05,
  extend = FALSE,
  n_boot = 2000,
  tau = NULL,
  weights = NULL
) {
  
  # Convert time and status to matrices.
  statuses <- as.matrix(statuses)
  times <- as.matrix(times)
  n_comp <- ncol(times)
  
  InputsChecks(
    statuses = statuses, 
    times = times
  )
  
  # Weights.
  if (is.null(weights)) {
    weights <- rep(1, n_comp)
  }
  
  # Truncation time.
  if (is.null(tau)) {
    tau <- MaxEventTime(statuses = statuses, times = times)
  }

  # Calculate test statistic.
  aucs <- lapply(seq_len(n_comp), function(i) {
    auc <- RMST(
      status = statuses[, i], 
      time = times[, i],
      extend = extend,
      tau = tau
    )
    return(auc)
  })
  aucs <- do.call(c, aucs)
  obs_stat <- sum(weights * aucs)
  
  # Calculate influence functions.
  psi <- lapply(seq_len(n_comp), function(i) {
    psi_i <- CalcPsiRMST(
      status = statuses[, i],
      time = times[, i],
      tau = tau
    )
    return(psi_i)
  })
  psi <- do.call(cbind, psi)
  psi <- psi %*% weights
  
  # Run perturbations. 
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

