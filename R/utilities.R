# Purpose: Utility functions.

#' Maximum Event Time
#' 
#' Find the maximum time with a status == 1.
#' 
#' @param statuses Array with 1 or more status variables.
#' @param times Array with 1 or more event-time variables.
#' @return Numeric maximum event time.
#' @export
MaxEventTime <- function(
  statuses,
  times
) {
  
  # Convert time and status to matrices.
  statuses <- as.matrix(statuses)
  times <- as.matrix(times)
  n_times <- ncol(times)
  
  max_event_times <- lapply(seq_len(n_times), function(i) {
    current_status <- statuses[, i]
    current_time <- times[, i]
    return(max(current_time[current_status == 1]))
  })
  max_event_times <- do.call(c, max_event_times)
  max_event_times <- max(0, max_event_times)
  return(max(max_event_times))
}

