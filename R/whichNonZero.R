#' Find non-zero entries of a matrix
#'
#' Finds the non-zero entries of a matrix in the most efficient manner for each matrix representation.
#' Not sure there's much more to say here.
#' 
#' @param x A numeric matrix-like object, usually sparse in content if not in representation.
#' @param BPPARAM A BiocParallelParam object from the \pkg{BiocParallel} package controlling how parallelization should be performed.
#' Only used when \code{x} is a \linkS4class{DelayedMatrix} object; defaults to no parallelization.
#' @param ... For the generic, additional arguments to pass to the specific methods.
#'
#' For the methods, additional arguments that are currently ignored.
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
#' @import methods
setGeneric("whichNonZero", function(x, ...) standardGeneric("whichNonZero"))

#' @export
#' @rdname whichNonZero
#' @importFrom DelayedArray which
setMethod("whichNonZero", "ANY", function(x, ...) {
    idx <- which(x!=0, arr.ind=TRUE)
    list(i=idx[,1], j=idx[,2], x=x[idx])
})

#' @export
#' @rdname whichNonZero
#' @importClassesFrom Matrix TsparseMatrix
setMethod("whichNonZero", "TsparseMatrix", function(x, ...) {
    list(i=x@i+1L, j=x@j+1L, x=x@x)
})

#' @export
#' @rdname whichNonZero
#' @importClassesFrom Matrix CsparseMatrix
setMethod("whichNonZero", "CsparseMatrix", function(x, ...) {
    i <- x@i + 1L
    d <- diff(x@p)
    j <- rep.int(seq_along(d), d)
    list(i=i, j=j, x=x@x)
})

#' @export
#' @rdname whichNonZero
#' @importClassesFrom DelayedArray SparseArraySeed
#' @importFrom DelayedArray nzindex nzdata
setMethod("whichNonZero", "SparseArraySeed", function(x, ...) {
    idx <- nzindex(x)
    if (ncol(idx)!=2) {
        stop("'x' should be a 2-dimensional SparseArraySeed")
    }
    list(i=idx[,1], j=idx[,2], x=nzdata(x))
})

#' @export
#' @rdname whichNonZero
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM 
#' @importClassesFrom DelayedArray SparseArraySeed
setMethod("whichNonZero", "DelayedMatrix", function(x, BPPARAM=NULL, ...) {
    oldBP <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldBP))

    out <- as(x, "SparseArraySeed")
    callGeneric(out)
})
