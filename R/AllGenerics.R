#' Initialize matrix in C++ memory space
#'
#' Initialize a \pkg{tatami} matrix object in C++ memory space from an abstract numeric R matrix.
#' This object simply references the R memory space and avoids making any copies of its own, so it can be cheaply re-created when needed inside each function.
#'
#' @param x A matrix-like object, typically from the \pkg{Matrix} or \pkg{DelayedArray} packages.
#' Alternatively, an external pointer from a previous call to \code{initializeCpp}, which is returned without modification.
#' @param ... Further arguments used by specific methods, such as:
#' \itemize{
#' \item \code{.check.na}, a logical vector indicating whether to check for \code{NA} values in integer and logical matrices.
#' If \code{TRUE} (the default), any \code{NA}s are cast to their double-precision equivalents when reading from the tatami matrix.
#' This can be set to \code{FALSE} to improve performance if the caller does not care about correctly handling \code{NA}s in \code{x}.
#' }
#' Fields should generally be prefixed by the matrix type, to avoid conflicts with arguments from other packages.
#' For example, \code{hdf5.realize} can be used in \pkg{beachmat.hdf5} to load a HDF5-backed matrix into memory.
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
#' initializeCpp
#' initializeCpp,ANY-method
#' initializeCpp,matrix-method
#' initializeCpp,externalptr-method
#' initializeCpp,dgeMatrix-method
#' initializeCpp,lgeMatrix-method
#' initializeCpp,dgCMatrix-method
#' initializeCpp,dgRMatrix-method
#' initializeCpp,lgCMatrix-method
#' initializeCpp,lgRMatrix-method
#' initializeCpp,SVT_SparseMatrix-method
#' initializeCpp,DelayedMatrix-method
#' initializeCpp,DelayedAbind-method
#' initializeCpp,DelayedAperm-method
#' initializeCpp,DelayedSubset-method
#' initializeCpp,DelayedSetDimnames-method
#' initializeCpp,DelayedUnaryIsoOpWithArgs-method
#' initializeCpp,DelayedUnaryIsoOpStack-method
#' initializeCpp,DelayedNaryIsoOp-method
#' @import methods
setGeneric("initializeCpp", function(x, ...) standardGeneric("initializeCpp"))
