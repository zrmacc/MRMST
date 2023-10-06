# Purpose: Two-sample multiple RMST calculation.
# Updated: 2023-10-06


#' Partition Data
#' 
#' @param arm Treatment arm, coded as 0 for reference, 1 for treatment.
#' @param statuses Array with 1 or more status variables.
#' @param times Array with 1 or more event-time variables.
#' @return List of partitioned data.
#' @noRd 
PartitionData <- function(arm, statuses, times) {
  
  statuses0 <- statuses[arm == 0, , drop = FALSE]
  statuses1 <- statuses[arm == 1, , drop = FALSE]
  
  times0 <- times[arm == 0, , drop = FALSE]
  times1 <- times[arm == 1, , drop = FALSE]
  
  out <- list(
    statuses0 = statuses0,
    statuses1 = statuses1,
    times0 = times0,
    times1 = times1
  )
  return(out)
}


#' Add Arm
#' 
#' Add treatment arm to OneSample results table.
#' 
#' @param arm Treatment arm.
#' @param one_sample OneSample object.
#' @return OneSample object.
#' @noRd
AddArm <- function(arm, one_sample) {
  tab <- one_sample@Table
  tab$arm <- arm
  tab <- tab %>% dplyr::relocate("arm")
  one_sample@Table <- tab
  return(one_sample)
}



#' Difference
#'
#' @param one_sample Data.frame of one sample results.
#' @param alpha Type I error.
#' @return Data.frame for the difference.
#' @noRd
Diff <- function(one_sample, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  arm <- auc <- est <- se <- NULL
  out <- one_sample %>%
    dplyr::summarise(
      stat = "A1-A0",
      est = auc[arm == 1] - auc[arm == 0],
      se = sqrt(se[arm == 1]^2 + se[arm == 0]^2)
    ) %>%
    dplyr::mutate(
      lower = est - z * se,
      upper = est + z * se,
      p = 2 * stats::pnorm(
        q = abs(est) / se,
        lower.tail = FALSE
      )
    )
  return(out)
}


#' Ratio
#'
#' @param one_sample Data.frame of one sample results.
#' @param alpha Type I error.
#' @return Data.frame for the difference.
#' @noRd
Ratio <- function(one_sample, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  arm <- auc <- est <- log_se <- se <- NULL
  out <- one_sample %>%
    dplyr::summarise(
      stat = "A1/A0",
      est = auc[arm == 1] / auc[arm == 0],
      log_se = sqrt(se[arm == 1]^2 / auc[arm == 1]^2 + se[arm == 0]^2 / auc[arm == 0]^2)
    ) %>% 
    dplyr::mutate(
      lower = est * exp(-z * log_se),
      upper = est * exp(+z * log_se),
      se = est * log_se,
      p = 2 * stats::pnorm(
        q = abs(log(est)) / log_se,
        lower.tail = FALSE
      )
    ) %>% 
    dplyr::select(-log_se)
  return(out)
}


#' Difference and Ratio
#' 
#' @param arm0 Data.frame of arm 0 results.
#' @param arm1 Data.frame of arm 1 results.
#' @param alpha Type I error.
#' @return Data.frame of contrasts.
#' @noRd
DiffRatio <- function(arm0, arm1, alpha = 0.05) {
  
  one_sample <- rbind(arm0, arm1)
  one_sample$arm <- c(0, 1)
  
  diff <- Diff(one_sample, alpha = alpha)
  ratio <- Ratio(one_sample, alpha = alpha)
  
  out <- rbind(diff, ratio)
  return(out)
}


#' Two Sample Calculation
#' 
#' @param arm Treatment arm, coded as 0 for reference, 1 for treatment.
#' @param statuses Array with 1 or more status variables.
#' @param times Array with 1 or more event-time variables.
#' @param alpha Type I error.
#' @param extend Extend the Kaplan-Meier curves if `tau` is set beyond the
#'   maximum observation time?
#' @param n_boot Number of perturbations for variance calculation.
#' @param tau Truncation time.
#' @param weights Optional weight vector. Defaults to equal weights.
#' @export
TwoSample <- function(
    arm,
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
  n_times <- ncol(times)

  InputsChecks(statuses = statuses, times = times, arm = arm)
  
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
  
  # Partition data.
  split_data <- PartitionData(arm = arm, statuses = statuses, times = times)
  
  # Truncation time.
  if (is.null(tau)) {
    tau0 <- MaxEventTime(statuses = split_data$statuses0, times = split_data$times0)
    tau1 <- MaxEventTime(statuses = split_data$statuses1, times = split_data$times1)
    tau <- min(tau0, tau1)
  }
  
  # Per arm results.
  arm0 <- OneSample(
    statuses = split_data$statuses0,
    times = split_data$times0,
    alpha = alpha,
    extend = extend,
    n_boot = n_boot,
    tau = tau,
    weights = weights
  )
  arm0 <- AddArm(arm = 0, one_sample = arm0)
  
  arm1 <- OneSample(
    statuses = split_data$statuses1,
    times = split_data$times1,
    alpha = alpha,
    extend = extend,
    n_boot = n_boot,
    tau = tau,
    weights = weights
  )
  arm1 <- AddArm(arm = 1, one_sample = arm1)
  
  # Contrast treatment arms.
  contrast <- DiffRatio(
    arm0 = arm0@Table,
    arm1 = arm1@Table,
    alpha = alpha
  )
  
  # Output.
  out <- methods::new(
    "TwoSample",
    Arm0 = arm0,
    Arm1 = arm1,
    Contrast = contrast
  )
  return(out)
}


