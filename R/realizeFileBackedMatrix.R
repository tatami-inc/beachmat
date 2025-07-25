#' Realize a file-backed DelayedMatrix
#'
#' Realize a file-backed DelayedMatrix into its corresponding in-memory format.
#'
#' @param x A \link[DelayedArray]{DelayedMatrix} object.
#'
#' @return 
#' For \code{realizeFileBackedMatrix}, an ordinary matrix or a dgCMatrix, depending on whether \code{\link[DelayedArray]{is_sparse}(x)}.
#'
#' For \code{isFileBackedMatrix}, a logical scalar indicating whether \code{x} has file-backed components.
#'
#' @details
#' A file-backed matrix representation is recognized based on whether it has a \code{\link[DelayedArray]{path}} method for any one of its seeds. 
#' If so, and the \code{"beachmat.realizeFileBackedMatrix"} option is not \code{FALSE}, we will load it into memory.
#' This is intended for DelayedMatrix objects that have already been subsetted (e.g., to highly variable genes),
#' which can be feasibly loaded into memory for rapid calculations.
#'
#' @author Aaron Lun
#'
#' @examples
#' mat <- matrix(rnorm(50), ncol=5)
#' realizeFileBackedMatrix(mat) # no effect
#' 
#' library(HDF5Array)
#' mat2 <- as(mat, "HDF5Array")
#' realizeFileBackedMatrix(mat2) # realized into memory
#' 
#' @export
#' @importFrom DelayedArray is_sparse
#' @importClassesFrom Matrix sparseMatrix
realizeFileBackedMatrix <- function(x) {
    if (isFileBackedMatrix(x) && getOption("beachmat.realizeFileBackedMatrix", TRUE)) {
        if (is_sparse(x)) {
            x <- as(x, "sparseMatrix")
        } else {
            x <- as.matrix(x)
        }
    }

    x    
}

#' @export
#' @rdname realizeFileBackedMatrix
#' @importFrom BiocGenerics path
#' @importFrom DelayedArray seedApply 
#' @importClassesFrom DelayedArray DelayedMatrix
isFileBackedMatrix <- function(x) {
    if (!is(x, "DelayedMatrix")) {
        return(FALSE)
    }

    # Figure out if any of the underlying seeds have a path() method.
    has.path <- tryCatch({
        seedApply(x, function(i) {
            p <- try(path(i), silent=TRUE)
            !is(p, "try-error")
        })
    }, error=function(e) {
        FALSE
    })

    any(unlist(has.path))
}
