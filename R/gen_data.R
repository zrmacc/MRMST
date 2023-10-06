# Purpose: Generate data.
# Updated: 2023-09-26


#' Generate Linear Predictor
#' 
#' @param beta_event Numeric vector of log rate ratios for the event rate.
#' @param covariates Numeric design matrix.
#' @param n_subj Number of subjects.
#' @return Numeric vector.
#' @noRd
GenLinPred <- function(beta_event, covariates, n_subj) {
  if (!is.null(beta_event) & !is.null(covariates)) {
    covariates <- data.matrix(covariates)
    eta <- as.numeric(covariates %*% beta_event)
  } else {
    eta <- rep(0, n_subj)
  }
  return(eta)
}


#' Generate Frailty
#'
#' @param frailty_variance Variance of the gamma frailty.
#' @param n_subj Number of subjects.
#' @return Numeric vector
#' @noRd
GenFrailty <- function(frailty_variance, n_subj) {
  if (frailty_variance > 0) {
    theta <- 1 / frailty_variance
    frailty <- stats::rgamma(n = n_subj, shape = theta, rate = theta)
  } else {
    frailty <- rep(1, n_subj)
  }
  return(frailty)
}


#' Generate Event Time Data
#' 
#' Generates time and status for `n_event` parallel events.
#' 
#' @param df Data frame including the 
#' @param n_event Number of parallel events.
#' @param tau Truncation time.
#' @return Data frame.
#' @noRd
GenEventTimeData <- function(df, n_event, tau) {
  
  event_rate <- df$event_rate
  censor_rate <- df$censor_rate
  n_subj <- nrow(df)
  
  # Generate event times.
  event_time <- lapply(seq_len(n_event), function(i) {
    event_time <- stats::rexp(n = n_subj, rate = event_rate)
    censor_time <- stats::rexp(n = n_subj, rate = censor_rate)
    time <- pmin(event_time, censor_time, tau)
    status <- 1 * (time == event_time)
    out <- data.frame(time = time, status = status)
    return(out)
  })
  event_time <- do.call(cbind, event_time)
  
  # Apply column names.
  col_names <- lapply(seq_len(n_event), function(i) {
    return(paste0(c("time", "status"), i))
  })
  col_names <- do.call(c, col_names)
  colnames(event_time) <- col_names
  
  return(event_time)
}


#' Generate Data
#' 
#' Generate event-time data for multiple parallel events, possibly depending on
#' covariates.
#' 
#' @param base_event_rate Baseline arrival rate for events.
#' @param beta_event Numeric vector of log rate ratios for the event rate.
#' @param censor_rate Arrival rate for the censoring time.
#' @param covariates Numeric design matrix.
#' @param frailty_variance Variance of the gamma frailty.
#' @param min_event_rate Minimum subject-specific event rate. Must be positive.
#' @param n_event Number of parallel events.
#' @param n_subj Number of subjects. Overwritten by the number of rows of the
#'   of the covariate matrix, if covariates are provided.
#' @param tau Truncation time.
#' @return Data.frame.
#' @export
GenData <- function(
  base_event_rate = 1.0,
  beta_event = NULL,
  censor_rate = 0.25,
  covariates = NULL,
  frailty_variance = 0.0,
  min_event_rate = 0.05,
  n_event = 2,
  n_subj = 100,
  tau = 4.0
) {
  
  # Linear predictor.
  eta <- GenLinPred(beta_event = beta_event, covariates = covariates, n_subj = n_subj)
  n_subj <- length(eta)
  df <- data.frame(idx = seq_len(n_subj))
  
  # Add covaraites, if provided.
  if (!is.null(covariates)) {
    df <- cbind(df, covariates)
  }
  
  # Calculate subject-specific event rate.
  df$event_rate <- base_event_rate * exp(eta)
  df$event_rate <- pmax(df$event_rate, min_event_rate)
  
  # Apply frailty.
  df$frailty <- GenFrailty(frailty_variance = frailty_variance, n_subj = n_subj)
  df$event_rate <- df$frailty * df$event_rate 
  
  # Add censoring rate.
  df$censor_rate <- censor_rate
  
  # Generate event times.
  event_times <- GenEventTimeData(df = df, n_event = n_event, tau = tau)
  out <- cbind(df, event_times)
  return(out)
}
