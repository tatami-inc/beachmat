#' Apply over blocks of columns or rows
#'
#' Apply a function over blocks of columns or rows using \pkg{DelayedArray}'s block processing mechanism.
#'
#' @param x A matrix-like object to be split into blocks and looped over.
#' This can be of any class that respects the matrix contract.
#' @param FUN A function that operates on columns or rows in \code{x},
#' for \code{colBlockApply} and \code{rowBlockApply} respectively.
#' Ordinary matrices, *gCMatrix or \linkS4class{SparseArraySeed} objects may be passed as the first argument.
#' @param ... Further arguments to pass to \code{FUN}.
#' @param grid An \linkS4class{ArrayGrid} object specifying how \code{x} should be split into blocks.
#' For \code{colBlockApply} and \code{rowBlockApply}, blocks should consist of consecutive columns and rows, respectively.
#' @param BPPARAM A BiocParallelParam object from the \pkg{BiocParallel} package,
#' specifying how parallelization should be performed across blocks.
#'
#' @return
#' A list of length equal to the number of blocks, 
#' where each entry is the output of \code{FUN} for the results of processing each the rows/columns in the corresponding block.
#'
#' @details
#' This is a wrapper around \code{\link{blockApply}} that is dedicated to looping across rows or columns of \code{x}.
#' The aim is to provide a simpler interface for the common task of \code{\link{apply}}ing across a matrix,
#' along with a few modifications to improve efficiency for parallel processing and for natively supported \code{x}.
#' Note that the fragmentation of \code{x} into blocks is not easily predictable, 
#' meaning that \code{FUN} should be capable of operating on each row/column independently.
#' 
#' @examples
#' x <- matrix(runif(10000), ncol=10)
#' str(colBlockApply(x, colSums))
#' str(rowBlockApply(x, rowSums))
#' 
#' library(Matrix)
#' y <- rsparsematrix(10000, 10000, density=0.01)
#' str(colBlockApply(y, colSums))
#' str(rowBlockApply(y, rowSums))
#' 
#' library(DelayedArray)
#' z <- DelayedArray(y) + 1
#' str(colBlockApply(z, colSums))
#' str(rowBlockApply(z, rowSums))
#'
#' # We can also force multiple blocks:
#' library(BiocParallel)
#' BPPARAM <- SnowParam(2)
#' str(colBlockApply(x, colSums, BPPARAM=BPPARAM))
#' str(rowBlockApply(x, rowSums, BPPARAM=BPPARAM))
#'
#' @seealso
#' \code{\link{blockApply}}, for the original \pkg{DelayedArray} implementation.
#' 
#' @export
#' @importFrom DelayedArray getAutoBPPARAM
colBlockApply <- function(x, FUN, ..., grid=NULL, BPPARAM=getAutoBPPARAM()) {
    .blockApply2(x, FUN=FUN, ..., grid=grid, BPPARAM=BPPARAM, by.row=FALSE)
}

#' @export
#' @rdname colBlockApply
#' @importFrom DelayedArray getAutoBPPARAM
rowBlockApply <- function(x, FUN, ..., grid=NULL, BPPARAM=getAutoBPPARAM()) {
    .blockApply2(x, FUN=FUN, ..., grid=grid, BPPARAM=BPPARAM, by.row=TRUE)
}

#' @importFrom methods is
#' @importFrom DelayedArray blockApply rowAutoGrid colAutoGrid
#' makeNindexFromArrayViewport getAutoBlockLength type
.blockApply2 <- function(x, FUN, ..., grid, BPPARAM, by.row=FALSE) {
    native <- is.matrix(x) || is(x, "lgCMatrix") || is(x, "dgCMatrix") 

    if (is.null(grid)) {
        nworkers <- if (is.null(BPPARAM)) 1L else BiocParallel::bpnworkers(BPPARAM)

        if (!native || nworkers != 1L) {
            # Scaling down the block length so that each worker is more likely to get a task.
            max.block.length <- getAutoBlockLength(type(x))
            expected.block.length <- length(x) / nworkers
            block.length <- min(max.block.length, expected.block.length)

            if (by.row) {
                grid <- rowAutoGrid(x, block.length=block.length)
            } else {
                grid <- colAutoGrid(x, block.length=block.length)
            }
        }
    }

    if (!native) {
        return(blockApply(x, FUN=FUN, ..., grid=grid, as.sparse=NA, BPPARAM=BPPARAM))

    } else if (is.null(grid) || length(grid)==1L) {
        # Avoid overhead of block processing if there isn't any grid.
        if (is.null(grid)) {
            grid <- RegularArrayGrid(dim(x))
        }
        attr(x, "from_grid") <- grid
        attr(x, "block_id") <- 1L
        return(list(FUN(x, ...)))

    }  else {
        # Break up the native matrix in the parent to ensure that we only 
        # need to serialize the chunks to the child. Note that 'fragments' 
        # still contains objects in their native format.
        fragments <- vector("list", length(grid))
        common.args <- list(x=x, drop=FALSE)

        for (i in seq_along(fragments)) {
            idx <- makeNindexFromArrayViewport(grid[[i]], expand.RangeNSBS=TRUE)

            if (is.null(idx[[1]])) {
                idx[[1]] <- substitute()
            }
            if (is.null(idx[[2]])) {
                idx[[2]] <- substitute()
            }
            names(idx) <- c("i", "j")

            block <- do.call("[", c(common.args, idx))
            attr(block, "from_grid") <- grid
            attr(block, "block_id") <- i
            fragments[[i]] <- block
        }

        DelayedArray:::bplapply2(fragments, FUN, ..., BPPARAM = BPPARAM)
    }
}
