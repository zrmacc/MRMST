# Purpose: Input checks.
# Updated: 2023-10-06

#' MRMST Input Checks
#' 
#' @param statuses Array with 1 or more status variables.
#' @param times Array with 1 or more event-time variables.
#' @param arm Treatment arm.
#' @noRd
InputsChecks <- function(statuses, times, arm = NULL) {
  
  n_statuses <- ncol(statuses)
  n_times <- ncol(times)
  if (n_statuses != n_times) {
    stop("Number of status variables does not match the number of time variables.")
  }
  
  if (!is.null(arm)) {
    arm_levels <- sort(unique(arm))
    if (!all(arm_levels == c(0, 1))) {
      stop("Arm should have exactly two levels, 0 for reference, 1 for treatment.")
    }
  }
  
}

