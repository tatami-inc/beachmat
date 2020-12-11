# This checks the sparse row subsetting functions.
# library(testthat); library(beachmat); source("test-sparse.R")

chunk_by_row_fast <- function(x, grid) {
    beachmat:::.prepare_sparse_row_subset(x, grid)
    lapply(grid, FUN=beachmat:::.subset_matrix, x=x)
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
