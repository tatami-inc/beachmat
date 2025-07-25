#' Apply over blocks of columns or rows
#'
#' Apply a function over blocks of columns or rows using \pkg{DelayedArray}'s block processing mechanism.
#'
#' @param x A matrix-like object to be split into blocks and looped over.
#' This can be of any class that respects the matrix contract.
#' @param FUN A function that operates on columns or rows in \code{x},
#' for \code{colBlockApply} and \code{rowBlockApply} respectively.
#' Ordinary matrices, CsparseMatrix or \link[SparseArray]{SparseMatrix} objects may be passed as the first argument.
#' @param ... Further arguments to pass to \code{FUN}.
#' @param grid An \link[DelayedArray]{ArrayGrid} object specifying how \code{x} should be split into blocks.
#' For \code{colBlockApply} and \code{rowBlockApply}, blocks should consist of consecutive columns and rows, respectively.
#' Alternatively, this can be set to \code{TRUE} or \code{FALSE}, see Details.
#' @param coerce.sparse Logical scalar indicating whether blocks of a sparse \link[DelayedArray]{DelayedMatrix} \code{x} should be automatically coerced into CsparseMatrix objects. 
#' @param BPPARAM A BiocParallelParam object from the \pkg{BiocParallel} package,
#' specifying how parallelization should be performed across blocks.
#'
#' @return
#' A list of length equal to the number of blocks, 
#' where each entry is the output of \code{FUN} for the results of processing each the rows/columns in the corresponding block.
#'
#' @details
#' This is a wrapper around \code{\link[DelayedArray]{blockApply}} that is dedicated to looping across rows or columns of \code{x}.
#' The aim is to provide a simpler interface for the common task of \code{\link{apply}}ing across a matrix,
#' along with a few modifications to improve efficiency for parallel processing and for natively supported \code{x}.
#'
#' Note that the fragmentation of \code{x} into blocks is not easily predictable,
#' meaning that \code{FUN} should be capable of operating on each row/column independently.
#' Users can retrieve the current location of each block of \code{x} by calling \code{\link[DelayedArray]{currentViewport}} inside \code{FUN}.
#'
#' If \code{grid} is not explicitly set to an \link[DelayedArray]{ArrayGrid} object, it can take several values:
#' \itemize{
#' \item If \code{TRUE}, the function will choose a grid that (i) respects the memory limits in \code{\link[DelayedArray]{getAutoBlockSize}}
#' and (ii) fragments \code{x} into sufficiently fine chunks that every worker in \code{BPPARAM} gets to do something.
#' If \code{FUN} might make large allocations, this mode should be used to constrain memory usage.
#' \item The default \code{grid=NULL} is very similar to \code{TRUE} 
#' except that that memory limits are ignored when \code{x} is of any type that can be passed directly to \code{FUN}.
#' This avoids unnecessary copies of \code{x} and is best used when \code{FUN} itself does not make large allocations.
#' \item If \code{FALSE}, the function will choose a grid that covers the entire \code{x}.
#' This is provided for completeness and is only really useful for debugging.
#' }
#'
#' The default of \code{coerce.sparse=TRUE} will generate dgCMatrix objects during block processing of a sparse DelayedMatrix \code{x}.
#' This is convenient as it avoids the need for \code{FUN} to specially handle \link[SparseArray]{SparseMatrix} objects from the \pkg{SparseArray} package.
#' If the coercion is not desired (e.g., to preserve integer values in \code{x}), it can be disabled with \code{coerce.sparse=FALSE}.
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
#' \code{\link[DelayedArray]{blockApply}}, for the original \pkg{DelayedArray} implementation.
#'
#' \code{\link{toCsparse}}, to convert SparseMatrix objects to CsparseMatrix objects prior to further processing in \code{FUN}.
#' 
#' @export
#' @importFrom DelayedArray getAutoBPPARAM
colBlockApply <- function(x, FUN, ..., grid=NULL, coerce.sparse=TRUE, BPPARAM=getAutoBPPARAM()) {
    .blockApply2(x, FUN=FUN, ..., grid=grid, coerce.sparse=coerce.sparse, BPPARAM=BPPARAM, beachmat_by_row=FALSE)
}

#' @export
#' @rdname colBlockApply
#' @importFrom DelayedArray getAutoBPPARAM
rowBlockApply <- function(x, FUN, ..., grid=NULL, coerce.sparse=TRUE, BPPARAM=getAutoBPPARAM()) {
    .blockApply2(x, FUN=FUN, ..., grid=grid, coerce.sparse=coerce.sparse, BPPARAM=BPPARAM, beachmat_by_row=TRUE)
}

#' @importFrom methods is
#' @importFrom DelayedArray blockApply DummyArrayGrid isPristine seed DelayedArray
.blockApply2 <- function(x, FUN, ..., grid, BPPARAM, coerce.sparse=TRUE, beachmat_by_row=FALSE) {
    if (is(x, "DelayedArray")) {
        if (isPristine(x)) {
            cur.seed <- seed(x)
            output <- .sanitize_if_native(cur.seed)
        } else {
            output <- list(native=FALSE, csparse=FALSE, x=x)
        }
    } else {
        output <- .sanitize_if_native(x)
    }

    native <- output$native
    csparse <- output$csparse
    x <- output$x

    nworkers <- if (is.null(BPPARAM)) 1L else BiocParallel::bpnworkers(BPPARAM)

    if (isFALSE(grid)) {
        grid <- DummyArrayGrid(dim(x))
    }
    
    if (!native) {
        if (!is(grid, "ArrayGrid")) {
            grid <- .define_multiworker_grid(x, nworkers, beachmat_by_row=beachmat_by_row) 
        }
        if (coerce.sparse) {
            output <- blockApply(x, FUN=.sparse_helper, beachmat_internal_FUN=FUN, ..., grid=grid, as.sparse=NA, BPPARAM=BPPARAM)
        } else {
            output <- blockApply(x, FUN=FUN, ..., grid=grid, as.sparse=NA, BPPARAM=BPPARAM)
        }

    } else {
        if (!is(grid, "ArrayGrid")) {
            if (isTRUE(grid)) {
                grid <- .define_multiworker_grid(x, nworkers, beachmat_by_row=beachmat_by_row)
            } else if (nworkers!=1) {
                grid <- .define_multiworker_grid(x, nworkers, beachmat_by_row=beachmat_by_row, 
                    max.block.length=.Machine$integer.max)
            } else {
                grid <- DummyArrayGrid(dim(x))
            }
        }

        if (length(grid)==1L) {
            # Avoid overhead of block subsetting if there isn't any grid.
            frag.info <- list(grid, 1L, x)
            output <- list(.viewport_helper(frag.info, beachmat_internal_FUN=FUN, ...))

        } else {
            if (beachmat_by_row && csparse) {
                # This predefines the indices for faster chunking later on:
                # try rowBlockApply(y, rowSums, grid=TRUE) with and without this line.
                extra <- .prepare_sparse_row_subset(x, grid)
            } else {
                extra <- NULL
            }

            if (is.null(BPPARAM) || is(BPPARAM, "SerialParam") || is(BPPARAM, "MulticoreParam")) {
                # In serial or shared-memory cases, we can do the subsetting in each worker.
                # This avoids the effective copy of the entire matrix when we split it up,
                # while also bypassing any need to serialize the entire matrix to the workers.
                output <- DelayedArray:::bplapply2(seq_along(grid), function(i) {
                    block <- .subset_matrix(x, grid[[i]], extra[[i]])
                    frag.info <- list(grid, i, block)
                    .viewport_helper(frag.info, beachmat_internal_FUN=FUN, ...)
                }, BPPARAM=BPPARAM)

            } else {
                # Break up the native matrix in the parent to ensure that we only 
                # need to serialize the chunks to the child. Note that 'fragments' 
                # still contains objects in their native format.
                fragments <- vector("list", length(grid))
                for (i in seq_along(fragments)) {
                    block <- .subset_matrix(x, grid[[i]], extra[[i]])
                    fragments[[i]] <- list(grid, i, block)
                }
                output <- DelayedArray:::bplapply2(fragments, FUN=.viewport_helper, beachmat_internal_FUN=FUN, ..., BPPARAM=BPPARAM)
            }
        }
    }

    output
}

#' @importFrom DelayedArray rowAutoGrid colAutoGrid getAutoBlockLength type
.define_multiworker_grid <- function(x, nworkers, beachmat_by_row, 
    max.block.length=getAutoBlockLength(type(x)))
{
    # Scaling down the block length so that each worker is more likely to get a task.
    if (beachmat_by_row) {
        expected.block.length <- max(1, ceiling(nrow(x) / nworkers) * ncol(x))
        block.length <- min(max.block.length, expected.block.length)
        rowAutoGrid(x, block.length=block.length)
    } else {
        expected.block.length <- max(1, ceiling(ncol(x) / nworkers) * nrow(x))
        block.length <- min(max.block.length, expected.block.length)
        colAutoGrid(x, block.length=block.length)
    }
}

.prepare_sparse_row_subset <- function(x, grid) {
    nrows <- dims(grid)[,1]
    limits <- cumsum(nrows)
    grid <- fragment_sparse_rows(x@i, x@p, limits)

    last <- 0L
    for (i in seq_along(nrows)) {
        grid[[i]][[3]] <- c(nrows[i], limits[i])
        choice <- last + seq_len(nrows[i])
        grid[[i]][4] <- list(rownames(x)[choice]) # possibly NULL.
        last <- last + nrows[i]
    }

    grid
}

#' @importClassesFrom Matrix lgCMatrix dgCMatrix
.sanitize_if_native <- function(x) {
    native <- FALSE
    csparse <- FALSE

    if (is.matrix(x)) {
        native <- TRUE
    } else if (is(x, "dgCMatrix")) {
        native <- TRUE
        csparse <- TRUE
        if (class(x)[1] != "dgCMatrix") {
            x <- as(x, "dgCMatrix") # Avoid difficulties with subclasses.
        }
    } else if (is(x, "lgCMatrix")) {
        native <- TRUE
        csparse <- TRUE
        if (class(x)[1] != "lgCMatrix") {
            x <- as(x, "lgCMatrix") # Avoid difficulties with subclasses.
        }
    } 

    list(native=native, csparse=csparse, x=x)
}

#' @useDynLib beachmat
#' @importFrom Rcpp sourceCpp
#' @importFrom DelayedArray makeNindexFromArrayViewport
#' @importFrom methods new
.subset_matrix <- function(x, vp, rowsparse=NULL) {
    if (is.null(rowsparse)) {
        idx <- makeNindexFromArrayViewport(vp, expand.RangeNSBS=TRUE)
        i <- idx[[1]]
        j <- idx[[2]]
        if (!is.null(i) && !is.null(j)) {
            x[i, j, drop=FALSE]
        } else if (!is.null(j)) {
            x[, j, drop=FALSE]
        } else if (!is.null(i)) {
            x[i, , drop=FALSE]
        } else {
            x
        }
    } else {
        # This had damn better be a sparse matrix, subsetted by row!
        idx <- sparse_subset_index(rowsparse[[1]], rowsparse[[2]])
        new(class(x), 
            x=x@x[idx], 
            i=x@i[idx] - (rowsparse[[3]][2] - rowsparse[[3]][1]), # adjusting row indices for the new matrix.
            p=rowsparse[[2]],
            Dim=c(rowsparse[[3]][1], ncol(x)),
            Dimnames=list(rowsparse[[4]], colnames(x)))
    }
}

#' @importFrom DelayedArray set_grid_context
.viewport_helper <- function(X, beachmat_internal_FUN, ...) {
    set_grid_context(X[[1]], X[[2]], X[[1]][[X[[2]]]])
    beachmat_internal_FUN(X[[3]], ...)
}

#' @importFrom DelayedArray set_grid_context effectiveGrid currentBlockId currentViewport
.sparse_helper <- function(beachmat_internal_x, beachmat_internal_FUN, ...) {
    # Re-setting the grid context, just in case we moved to a different process without the previous settings.
    set_grid_context(effectiveGrid(), currentBlockId(), currentViewport())
    beachmat_internal_FUN(toCsparse(beachmat_internal_x), ...)
}
