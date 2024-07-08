#' Find non-zero entries of a matrix
#'
#' This function is soft-deprecated; users are advised to use \code{\link{nzwhich}} and \code{\link{nzvals}} instead.
#' 
#' @param x A numeric matrix-like object, usually sparse in content if not in representation.
#' @param ... Further arguments, ignored.
#'
#' @return A list containing \code{i}, an integer vector of the row indices of all non-zero entries;
#' \code{j}, an integer vector of the column indices of all non-zero entries;
#' and \code{x}, a (usually atomic) vector of the values of the non-zero entries.
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- Matrix::rsparsematrix(1e6, 1e6, 0.000001)
#' out <- whichNonZero(x)
#' str(out)
#' 
#' @seealso
#' \code{\link{which}}, obviously.
#'
#' @export
#' @importFrom SparseArray nzwhich nzvals
whichNonZero <- function(x, ...) {
    idx <- nzwhich(x, arr.ind=TRUE)
    vals <- nzvals(x)
    list(i=idx[,1], j=idx[,2], x=vals)
}
