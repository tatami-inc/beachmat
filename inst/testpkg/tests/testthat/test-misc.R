# Assorted tests that don't fit anywhere else.
# library(testthat); source("test-misc.R")

test_that("Output mode choices are okay", {
    expect_identical(beachtest:::check_output_mode("simple", simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(beachtest:::check_output_mode("sparse", simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(beachtest:::check_output_mode("RLE", simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(beachtest:::check_output_mode("Psymm", simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(beachtest:::check_output_mode("sparse", simplify=FALSE, preserve.zero=TRUE), "sparse")
    expect_identical(beachtest:::check_output_mode("sparse", simplify=FALSE, preserve.zero=FALSE), "HDF5")
    expect_identical(beachtest:::check_output_mode("RLE", simplify=FALSE, preserve.zero=FALSE), "HDF5")
    expect_identical(beachtest:::check_output_mode("Psymm", simplify=FALSE, preserve.zero=FALSE), "HDF5")
    expect_identical(beachtest:::check_output_mode("HDF5", simplify=FALSE, preserve.zero=FALSE), "HDF5")
})

# Checking that construction of empty matrices are okay.

library(Matrix)
library(HDF5Array)
test_that("empty matrices are okay", {
    for (mode in c("matrix", "dgCMatrix", "HDF5Array")) {
        # No columns
        ref <- as(matrix(0, 10, 0), mode)
        out <- .Call(beachtest:::cxx_test_numeric_output, ref, 3L, integer(0))
  
        if (isS4(ref)) { 
            expect_s4_class(out[[1]], mode)
        } else {
            expect_true(is(out[[1]], "matrix"))
        }

        expect_identical(ncol(out[[1]]), 0L)
        expect_identical(ncol(out[[2]]), 0L)
        expect_identical(nrow(out[[1]]), 10L)
        expect_identical(nrow(out[[2]]), 10L)

        # No rows
        ref <- as(matrix(0, 0, 10), mode)
        out <- .Call(beachtest:::cxx_test_numeric_output, ref, 3L, integer(0))

        if (isS4(ref)) { 
            expect_s4_class(out[[1]], mode)
        } else {
            expect_true(is(out[[1]], "matrix"))
        }

        expect_identical(nrow(out[[1]]), 0L)
        expect_identical(nrow(out[[2]]), 0L)
        expect_identical(ncol(out[[1]]), 10L)
        expect_identical(ncol(out[[2]]), 10L)
    }
})

# Checking random column slices behave correctly.

set.seed(23456)
check_col_slices <- function(FUN, ...) { 
    A <- FUN(...)

    test.mat <- as.matrix(A)
    dimnames(test.mat) <- NULL
    slice.start <- sample(ncol(A), nrow(A), replace=TRUE)
    slice.end <- pmin(ncol(A), slice.start + sample(10, nrow(A), replace=TRUE))
    
    out <- .Call(beachtest:::cxx_test_sparse_numeric_slice, A, cbind(slice.start, slice.end))
    ref <- vector('list', nrow(A))
    for (x in seq_along(ref)) { 
        ref[[x]] <- as.vector(A[x,slice.start[x]:slice.end[x]])
    }
    expect_identical(out, ref)
}

library(Matrix)
test_that("Sparse numeric indexing with slices is okay", {
    check_col_slices(FUN=rsparsematrix, nrow=100, ncol=20, density=0.2)
    check_col_slices(FUN=rsparsematrix, nrow=100, ncol=20, density=0.1)
    check_col_slices(FUN=rsparsematrix, nrow=100, ncol=50, density=0.2)
})

# Repeating with RLE matrix.

library(DelayedArray)
rFUN <- function(nr=15, nc=10, density=0.2) {
    as(as.matrix(rsparsematrix(nr, nc, density)), "RleArray")
}

test_that("RLE numeric indexing with slices is okay", {
    check_col_slices(FUN=rFUN, nr=100, nc=20, density=0.2)
    check_col_slices(FUN=rFUN, nr=100, nc=20, density=0.1)
    check_col_slices(FUN=rFUN, nr=100, nc=50, density=0.2)
})

