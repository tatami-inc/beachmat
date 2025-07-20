#' Convert a SparseMatrix to a CsparseMatrix
#'
#' Exactly what it says in the title.
#'
#' @param x Any object produced by block processing with \code{\link{colBlockApply}} or \code{\link{rowBlockApply}}.
#' This can be a matrix, sparse matrix or a \link[SparseArray]{SparseMatrix} object. 
#'
#' @return \code{x} is returned unless it is a \link[SparseArray]{SparseMatrix} object,
#' in which case an appropriate CsparseMatrix object is returned instead.
#'
#' @details
#' This is intended for use inside functions to be passed to \code{\link{colBlockApply}} or \code{\link{rowBlockApply}}.
#' The idea is to pre-process blocks for user-defined functions that don't know how to deal with SparseMatrix objects,
#' which is often the case for R-defined functions that do not benefit from \pkg{beachmat}'s C++ abstraction.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(SparseArray)
#' out <- COO_SparseArray(c(10, 10),
#'     nzcoo=cbind(1:10, sample(10)),
#'     nzdata=runif(10))
#' toCsparse(out)
#'
#' @export
toCsparse <- function(x) {
    if (is(x, "SparseMatrix")) {
        x <- as(x, "CsparseMatrix")
    }
    x
}
