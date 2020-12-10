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
#' Alternatively, this can be set to \code{TRUE} or \code{FALSE}, see Details.
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
#'
#' Note that the fragmentation of \code{x} into blocks is not easily predictable,
#' meaning that \code{FUN} should be capable of operating on each row/column independently.
#' Users can retrieve the current location of each block within \code{x} with \code{\link{currentViewport}} inside \code{FUN}.
#'
#' If \code{grid} is not explicitly set to an \linkS4class{ArrayGrid} object, it can take several values:
#' \itemize{
#' \item If \code{TRUE}, the function will choose a grid that (i) respects the memory limits in \code{\link{getAutoBlockSize}}
#' and (ii) fragments \code{x} into sufficiently fine chunks that every worker in \code{BPPARAM} gets to do something.
#' If \code{FUN} might make large allocations, this mode should be used to constrain memory usage.
#' \item The default \code{grid=NULL} is very similar to \code{TRUE} 
#' except that that memory limits are ignored when \code{x} is of any type that can be passed directly to \code{FUN}.
#' This avoids unnecessary copies of \code{x} and is best used when \code{FUN} itself does not make large allocations.
#' \item If \code{FALSE}, the function will choose a grid that covers the entire \code{x}.
#' This is provided for completeness and is only really useful for debugging.
#' }
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
    .blockApply2(x, FUN=FUN, ..., grid=grid, BPPARAM=BPPARAM, beachmat_by_row=FALSE)
}

#' @export
#' @rdname colBlockApply
#' @importFrom DelayedArray getAutoBPPARAM
rowBlockApply <- function(x, FUN, ..., grid=NULL, BPPARAM=getAutoBPPARAM()) {
    .blockApply2(x, FUN=FUN, ..., grid=grid, BPPARAM=BPPARAM, beachmat_by_row=TRUE)
}

#' @importFrom methods is
#' @importFrom DelayedArray blockApply rowAutoGrid colAutoGrid
#' makeNindexFromArrayViewport getAutoBlockLength type DummyArrayGrid
#' isPristine seed DelayedArray
.blockApply2 <- function(x, FUN, ..., grid, BPPARAM, beachmat_by_row=FALSE) {
    if (is(x, "DelayedArray") && isPristine(x)) {
        cur.seed <- seed(x)
        if (.is_native(cur.seed)) {
            x <- cur.seed
        }
    } 
    native <- .is_native(x) 
    nworkers <- if (is.null(BPPARAM)) 1L else BiocParallel::bpnworkers(BPPARAM)

    if (is.logical(grid) || is.null(grid)) {
        if (isTRUE(grid) || !native || nworkers != 1L) {
            # Avoid block size limits for native matrices when grid=TRUE.
            if (native && !isTRUE(grid)) {
                max.block.length <- .Machine$integer.max
            } else {
                max.block.length <- getAutoBlockLength(type(x))
            }

            # Scaling down the block length so that each worker is more likely to get a task.
            if (beachmat_by_row) {
                expected.block.length <- max(1, ceiling(nrow(x) / nworkers) * ncol(x))
                block.length <- min(max.block.length, expected.block.length)
                grid <- rowAutoGrid(x, block.length=block.length)
            } else {
                expected.block.length <- max(1, ceiling(ncol(x) / nworkers) * nrow(x))
                block.length <- min(max.block.length, expected.block.length)
                grid <- colAutoGrid(x, block.length=block.length)
            }
        } else { 
            grid <- DummyArrayGrid(dim(x))
        } 
    }

    if (!native) {
        blockApply(x, FUN=FUN, ..., grid=grid, as.sparse=NA, BPPARAM=BPPARAM)

    } else if (length(grid)==1L) {
        # Avoid overhead of block subsetting if there isn't any grid.
        frag.info <- list(grid, 1L, x)
        list(.helper(frag.info, beachmat_internal_FUN=FUN, ...))

    } else {
        # Break up the native matrix in the parent to ensure that we only 
        # need to serialize the chunks to the child. Note that 'fragments' 
        # still contains objects in their native format.
        fragments <- vector("list", length(grid))
        common.args <- list(x=x, drop=FALSE)

        # However, if we're not parallelizing, then we just execute in 
        # the loop and avoid making an effective copy of the entire matrix.
        nopar <- nworkers==1L

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
            frag.info <- list(grid, i, block)

            if (nopar) {
                fragments[i] <- list(.helper(frag.info, beachmat_internal_FUN=FUN, ...))
            } else {
                fragments[[i]] <- frag.info
            }
        }

        if (nopar) {
            fragments
        } else {
            DelayedArray:::bplapply2(fragments, FUN=.helper, beachmat_internal_FUN=FUN, ..., BPPARAM=BPPARAM)
        }
    }
}

#' @importClassesFrom Matrix lgCMatrix dgCMatrix
.is_native <- function(x) {
    is.matrix(x) || is(x, "lgCMatrix") || is(x, "dgCMatrix")
}

#' @importFrom DelayedArray set_grid_context
.helper <- function(X, beachmat_internal_FUN, ...) {
    set_grid_context(X[[1]], X[[2]])
    beachmat_internal_FUN(X[[3]], ...)
}
