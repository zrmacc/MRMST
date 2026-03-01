# Purpose: Input checks.
# Updated: 2023-10-06

#' MRMST Input Checks
#'
#' @param statuses Array with 1 or more status variables.
#' @param times Array with 1 or more event-time variables.
#' @param arm Treatment arm.
#' @param tau Truncation time (optional).
#' @param weights Weight vector (optional).
#' @param n_comp Number of outcome components (optional, for weight length check).
#' @noRd
InputsChecks <- function(statuses, times, arm = NULL, tau = NULL, weights = NULL, n_comp = NULL) {
  n_statuses <- ncol(statuses)
  n_times <- ncol(times)
  if (n_statuses != n_times) {
    stop("Number of status variables does not match the number of time variables.")
  }

  if (nrow(times) == 0L) {
    stop("At least one observation is required.")
  }

  if (any(is.na(statuses)) || any(is.na(times))) {
    stop("Missing values are not allowed in statuses or times.")
  }

  if (!is.null(arm)) {
    if (any(is.na(arm))) {
      stop("Missing values are not allowed in arm.")
    }
    arm_levels <- sort(unique(arm))
    if (!identical(arm_levels, c(0, 1))) {
      stop("Arm should have exactly two levels, 0 for reference, 1 for treatment.")
    }
  }

  if (!is.null(tau)) {
    if (!is.finite(tau) || tau <= 0) {
      stop("tau must be a finite positive value.")
    }
  }

  if (!is.null(weights) && !is.null(n_comp)) {
    if (length(weights) != n_comp) {
      stop("Length of weights must equal the number of outcome components.")
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("Weights must be finite and non-negative.")
    }
  }

  return(invisible(NULL))
}

