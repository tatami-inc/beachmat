# Checks for proper realization of unknown matrices.
# library(testthat); library(beachmat); source("test-unknown.R")

set.seed(1000)
smallmat <- matrix(runif(10000), 20, 500)
test_that("unknown setup methods work correctly", {
    out <- beachmat:::setupUnknownMatrix(smallmat)
    expect_identical(out[[1]], dim(smallmat))
    expect_identical(out[[2]], dim(smallmat))

	old <- getOption("DelayedArray.block.size")
	options(DelayedArray.block.size=5*5*8)
    out <- beachmat:::setupUnknownMatrix(smallmat)
    expect_identical(out[[1]], dim(smallmat))
    expect_identical(out[[2]], c(5L, 5L))
	options(DelayedArray.block.size=old)
})

library(Matrix)
sparsemat <- rsparsematrix(100, 20, density=0.1)
test_that("realization by range works correctly", {
    real <- beachmat:::realizeByRange(sparsemat, c(0, 5), c(1, 2))
    expect_equal(real, as.matrix(sparsemat[1:5,2:3]))

    real <- beachmat:::realizeByRange(sparsemat, c(1, 6), c(2, 5))
    expect_equal(real, as.matrix(sparsemat[2:7,3:7]))
    
	real <- beachmat:::realizeByRange(sparsemat, c(1, 6), c(2, 5), transpose=TRUE)
    expect_equal(real, t(as.matrix(sparsemat[2:7,3:7])))

    real <- beachmat:::realizeByRange(sparsemat, c(3, 10), c(0, 20))
    expect_equal(real, as.matrix(sparsemat[4:13,,drop=FALSE]))

    real <- beachmat:::realizeByRange(sparsemat, c(3, 10), c(0, 0))
    expect_equal(real, as.matrix(sparsemat[4:13,integer(0),drop=FALSE]))

    real <- beachmat:::realizeByRange(sparsemat, c(0, 0), c(0, 20))
    expect_equal(real, as.matrix(sparsemat[integer(0),,drop=FALSE]))
})

test_that("realization by range vs index works correctly", {
    real <- beachmat:::realizeByRangeIndex(sparsemat, c(0, 5), 3:8)
    expect_equal(real, as.matrix(sparsemat[1:5,3:8]))

    real <- beachmat:::realizeByRangeIndex(sparsemat, c(0, 5), integer(0))
    expect_equal(real, as.matrix(sparsemat[1:5,integer(0),drop=FALSE]))

    real <- beachmat:::realizeByRangeIndex(sparsemat, c(0, 0), 2:5)
    expect_equal(real, as.matrix(sparsemat[integer(0),2:5,drop=FALSE]))

    real <- beachmat:::realizeByIndexRange(sparsemat, 1:5, c(2, 5))
    expect_equal(real, as.matrix(sparsemat[1:5,3:7]))

    real <- beachmat:::realizeByIndexRange(sparsemat, integer(0), c(2, 5))
    expect_equal(real, as.matrix(sparsemat[integer(0), 3:7, drop=FALSE]))

    real <- beachmat:::realizeByIndexRange(sparsemat, 3:7, c(0, 0))
    expect_equal(real, as.matrix(sparsemat[3:7,integer(0), drop=FALSE]))
})
