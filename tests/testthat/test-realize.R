# This checks the realizeFileBackedMatrix function.
# library(testthat); library(beachmat); source("test-realize.R")

test_that("realizeFileBackedMatrix is a no-op for in-memory objects", {
    mat <- matrix(rnorm(50), ncol=5)
    expect_identical(mat, realizeFileBackedMatrix(mat))

    library(DelayedArray)
    x <- DelayedArray(mat)
    expect_identical(x, realizeFileBackedMatrix(x))

    x1 <- x+1
    expect_identical(x1, realizeFileBackedMatrix(x1))
    expect_identical(BiocGenerics::cbind(x, x), realizeFileBackedMatrix(BiocGenerics::cbind(x, x)))

    y <- Matrix::rsparsematrix(50, 10, density=0.2)
    expect_identical(y, realizeFileBackedMatrix(y))

    y2 <- DelayedArray(y) + 10
    expect_identical(y2, realizeFileBackedMatrix(y2))
})

test_that("realizeFileBackedMatrix works for HDF5Array objects", {
    library(HDF5Array)
    mat <- matrix(rnorm(50), ncol=5)
    mat2 <- as(mat, "HDF5Array")

    expect_identical(mat, realizeFileBackedMatrix(mat2))

    expect_identical(cbind(mat, mat), realizeFileBackedMatrix(BiocGenerics::cbind(mat2, mat2)))

    expect_identical(cbind(mat, mat), realizeFileBackedMatrix(BiocGenerics::cbind(mat2, DelayedArray(mat))))

    expect_identical(mat + 1, realizeFileBackedMatrix(mat2 + 1))
})

test_that("realizeFileBackedMatrix works for sparse derivatives", {
    library(HDF5Array)
    mat <- rsparsematrix(50, 10, 0.2)
    mat2 <- writeHDF5Array(mat, as.sparse=TRUE)

    expect_identical(mat, realizeFileBackedMatrix(mat2))

    expect_identical(unname(as.matrix(mat) + 1), unname(realizeFileBackedMatrix(mat2 + 1)))
})
