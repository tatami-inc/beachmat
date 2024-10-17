#' Convert a SparseArraySeed to a CsparseMatrix
#'
#' Exactly what it says in the title.
#'
#' @param x Any object produced by block processing with \code{\link{colBlockApply}} or \code{\link{rowBlockApply}}.
#' This can be a matrix, sparse matrix or a two-dimensional \linkS4class{SparseArraySeed}.
#'
#' @return \code{x} is returned unless it is one of the \pkg{DelayedArray} sparse matrix classes,
#' in which case an appropriate \linkS4class{CsparseMatrix} object is returned instead.
#'
#' @details 
#' This is intended for use inside functions to be passed to \code{\link{colBlockApply}} or \code{\link{rowBlockApply}}.
#' The idea is to pre-process blocks for user-defined functions that don't know how to deal with SparseArraySeed objects,
#' which is often the case for R-defined functions that do not benefit from \pkg{beachmat}'s C++ abstraction.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(DelayedArray)
#' out <- SparseArraySeed(c(10, 10), 
#'     nzindex=cbind(1:10, sample(10)),
#'     nzdata=runif(10))
#' toCsparse(out)
#'
#' @export
toCsparse <- function(x) {
    if (is(x, "SparseArraySeed") || is(x, "SparseArray")) {
        x <- as(x, "CsparseMatrix")
    }
    x
}
