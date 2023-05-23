#' Initialize matrix in C++ memory space
#'
#' Initialize a \pkg{tatami} matrix object in C++ memory space from an abstract numeric R matrix.
#' This object simply references the R memory space and avoids making any copies of its own, so it can be cheaply re-created when needed inside each function.
#'
#' @param x A matrix-like object, typically from the \pkg{Matrix} or \pkg{DelayedArray} packages.
#' @param ... Further arguments used by specific methods.
#'
#' @return An external pointer to a C++ object containing a tatami matrix.
#'
#' @details
#' Do not attempt to serialize the return value; it contains a pointer to external memory, and will not be valid after a save/load cycle.
#' Users should not be exposed to the returned pointers; rather, developers should call \code{initialize} at the start to obtain a C++ object for further processing.
#' As mentioned before, this initialization process is very cheap so there is no downside from just recreating the object within each function body.
#'
#' @examples
#' # Mocking up a count matrix:
#' x <- Matrix::rsparsematrix(1000, 100, 0.1)
#' y <- round(abs(x))
#'
#' stuff <- initializeCpp(y)
#' stuff
#' 
#' @export
#' @aliases
#' initializeCpp,dgCMatrix-method
#' initializeCpp,dgRMatrix-method
#' @import methods
setGeneric("initializeCpp", function(x, ...) standardGeneric("initializeCpp"))
