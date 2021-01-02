# This checks the sparse row subsetting functions.
# library(testthat); library(beachmat); source("test-sparse.R")

chunk_by_row_fast <- function(x, grid) {
    grid <- beachmat:::.prepare_sparse_row_subset(x, grid)
    lapply(grid, FUN=beachmat:::.subset_matrix, x=x, vp=NULL)
}

chunk_by_row_ref <- function(x, grid) {
    lapply(seq_along(grid), FUN=function(i) beachmat:::.subset_matrix(x, grid[[i]]))
}

test_that("fast row chunking works correctly", {
    y <- Matrix::rsparsematrix(1000, 100, density=0.01)

    grid <- DelayedArray::RegularArrayGrid(dim(y))
    expect_identical(
        chunk_by_row_fast(y, grid),
        chunk_by_row_ref(y, grid)
    )

    grid <- DelayedArray::RegularArrayGrid(dim(y), spacings=c(10, ncol(y)))
    expect_identical(
        chunk_by_row_fast(y, grid),
        chunk_by_row_ref(y, grid)
    )

    grid <- DelayedArray::RegularArrayGrid(dim(y), spacings=c(nrow(y)/3, ncol(y)))
    expect_identical(
        chunk_by_row_fast(y, grid),
        chunk_by_row_ref(y, grid)
    )

    grid <- DelayedArray::RegularArrayGrid(dim(y), spacings=c(nrow(y)/10, ncol(y)))
    expect_identical(
        chunk_by_row_fast(y, grid),
        chunk_by_row_ref(y, grid)
    )

    # Handles special cases.
    expect_identical(
        chunk_by_row_fast(y[,0], grid),
        chunk_by_row_ref(y[,0], grid)
    )

    grid <- DelayedArray::RegularArrayGrid(c(0, ncol(y)))
    expect_identical(
        chunk_by_row_fast(y[0,], grid),
        chunk_by_row_ref(y[0,], grid)
    )
})

test_that("fast row names are passed along correctly", {
    y <- Matrix::rsparsematrix(1000, 100, density=0.01)
    rownames(y) <- sprintf("Y%i", seq_len(nrow(y)))

    grid <- DelayedArray::RegularArrayGrid(dim(y))
    out <- chunk_by_row_fast(y, grid)
    expect_identical(rownames(y), unlist(lapply(out, rownames)))

    grid <- DelayedArray::RegularArrayGrid(dim(y), spacings=c(10, ncol(y)))
    out <- chunk_by_row_fast(y, grid)
    expect_identical(rownames(y), unlist(lapply(out, rownames)))

    grid <- DelayedArray::RegularArrayGrid(dim(y), spacings=c(nrow(y)/3, ncol(y)))
    out <- chunk_by_row_fast(y, grid)
    expect_identical(rownames(y), unlist(lapply(out, rownames)))
})

library(DelayedArray)
test_that("grid viewport is correctly passed", {
    y <- Matrix::rsparsematrix(1000, 100, density=0.01)
    setAutoBlockSize(ncol(y) * 8 * 10)

    out <- rowBlockApply(y, function(x) currentViewport(), grid=TRUE)
    expect_identical(length(out), 100L)

    ref <- rowBlockApply(DelayedArray(y), function(x) currentViewport(), grid=TRUE)
    expect_identical(out, ref)

    setAutoBlockSize()
})
