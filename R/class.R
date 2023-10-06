# Purpose: Class for MRMST output.
# Updated: 2023-10-06


# -----------------------------------------------------------------------------
# One Sample
# -----------------------------------------------------------------------------

#' One Sample
#'
#' @slot AUC Area under the curve.
#' @slot K Number of outcomes.
#' @slot SE Standard error.
#' @slot Table Table.
#' @slot Tau Truncation time.
#' @name OneSample-class
#' @rdname OneSample-class
#' @exportClass OneSample
setClass(
  Class = "OneSample",
  representation = representation(
    AUC = "numeric",
    K = "integer",
    SE = "numeric",
    Table = "data.frame",
    Tau = "numeric"
  )
)


#' Print Method for OneSample Object
#'
#' Print method for objects of class \code{OneSample}.
#'
#' @param x An object of class \code{OneSample}.
#' @param ... Unused.
#' @export
print.OneSample <- function(x, ...) {
  
  cat(glue::glue("Multiple RMST with {x@K} components."))
  cat("\n")
  tab <- x@Table %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, digits = 3)))
  print(tab)
  
}


#' Show Method for OneSample Object
#'
#' @param object An object of class \code{OneSample}.
#' @rdname OneSample-method
#' @importFrom methods show
setMethod(
  f = "show",
  signature = c(object = "OneSample"),
  definition = function(object) {
    print.OneSample(x = object)
  }
)


# -----------------------------------------------------------------------------
# Two Sample
# -----------------------------------------------------------------------------

#' Two Sample
#'
#' @slot Arm0 Results for the reference arm.
#' @slot Arm1 Results for the treatment arm.
#' @slot Contrast Contrasts. 
#' @name TwoSample-class
#' @rdname TwoSample-class
#' @exportClass TwoSample
setClass(
  Class = "TwoSample",
  representation = representation(
    Arm0 = "OneSample",
    Arm1 = "OneSample",
    Contrast = "data.frame"
  )
)


#' Print Method for TwoSample Object
#'
#' Print method for objects of class \code{TwoSample}.
#'
#' @param x An object of class \code{TwoSample}.
#' @param ... Unused.
#' @export
print.TwoSample <- function(x, ...) {

  cat("Arm 0:\n")
  print.OneSample(x@Arm0)
  cat("\n")
  
  cat("Arm 1:\n")
  print.OneSample(x@Arm1)
  cat("\n")
  
  cat("Contrast:\n")
  contrast <- x@Contrast %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, digits = 3)))
  print(contrast)
  
}


#' Show Method for TwoSample Object
#'
#' @param object An object of class \code{TwoSample}.
#' @rdname TwoSample-method
#' @importFrom methods show
setMethod(
  f = "show",
  signature = c(object = "TwoSample"),
  definition = function(object) {
    print.TwoSample(x = object)
  }
)

