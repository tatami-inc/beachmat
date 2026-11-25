# library(testthat); library(beachmat); source("test-tatami-multiply.R")

set.seed(1000)
test_that("matrix-vector multiplication works as expected", {
    x1 <- Matrix::rsparsematrix(1000, 100, 0.1)
    ptr1 <- initializeCpp(x1)

    vec <- runif(ncol(x1))
    expected <- as.vector(x1 %*% vec)
    expect_equal(tatami.multiply(ptr1, vec, right=TRUE, num.threads=1), expected)
    expect_equal(tatami.multiply(ptr1, vec, right=TRUE, num.threads=2), expected)

    vec <- runif(nrow(x1))
    expected <- as.vector(rbind(vec) %*% x1)
    expect_equal(tatami.multiply(ptr1, vec, right=FALSE, num.threads=1), expected)
    expect_equal(tatami.multiply(ptr1, vec, right=FALSE, num.threads=2), expected)
})

library(DelayedArray)
test_that("matrix multiplication works as expected (normal)", {
    x1 <- Matrix::rsparsematrix(200, 500, 0.1)
    ptr1 <- initializeCpp(x1)
    mat <- matrix(runif(ncol(x1) * 3), ncol = 3)
    expected <-  as.matrix(x1 %*% mat)

    expect_equal(tatami.multiply(ptr1, mat, right=TRUE, num.threads=1), expected)
    expect_equal(tatami.multiply(ptr1, mat, right=TRUE, num.threads=2), expected)

    # Trying a dense LHS matrix.
    d1 <- as.matrix(x1)
    dptr1 <- initializeCpp(d1)
    expect_equal(tatami.multiply(dptr1, mat, right=TRUE, num.threads=1), expected)
    expect_equal(tatami.multiply(dptr1, mat, right=TRUE, num.threads=2), expected)

    # Trying row-major sparse LHS matrix.
    tptr1 <- initializeCpp(t(DelayedArray(t(x1))))
    expect_equal(tatami.multiply(tptr1, mat, right=TRUE, num.threads=1), expected)
    expect_equal(tatami.multiply(tptr1, mat, right=TRUE, num.threads=2), expected)

    # Trying row-major dense LHS matrix.
    tdptr1 <- initializeCpp(t(DelayedArray(t(d1))))
    expect_equal(tatami.multiply(tdptr1, mat, right=TRUE, num.threads=1), expected)
    expect_equal(tatami.multiply(tdptr1, mat, right=TRUE, num.threads=2), expected)

    # Trying two matrices.
    ptr2 <- initializeCpp(mat)
    expect_equal(tatami.multiply(ptr1, ptr2, right=TRUE, num.threads=1), expected)
    expect_equal(tatami.multiply(ptr1, ptr2, right=TRUE, num.threads=2), expected)
})

test_that("matrix multiplication works as expected (transposed)", {
    x1 <- Matrix::rsparsematrix(100, 1000, 0.1)
    ptr1 <- initializeCpp(x1)
    mat <- matrix(runif(nrow(x1) * 3), nrow = 3)
    expected <-  as.matrix(mat %*% x1)

    expect_equal(tatami.multiply(ptr1, mat, right=FALSE, num.threads=1), expected)
    expect_equal(tatami.multiply(ptr1, mat, right=FALSE, num.threads=2), expected)

    # Trying a dense LHS matrix.
    d1 <- as.matrix(x1)
    dptr1 <- initializeCpp(d1)
    expect_equal(tatami.multiply(dptr1, mat, right=FALSE, num.threads=1), expected)
    expect_equal(tatami.multiply(dptr1, mat, right=FALSE, num.threads=2), expected)

    # Trying row-major sparse LHS matrix.
    tptr1 <- initializeCpp(t(DelayedArray(t(x1))))
    expect_equal(tatami.multiply(tptr1, mat, right=FALSE, num.threads=1), expected)
    expect_equal(tatami.multiply(tptr1, mat, right=FALSE, num.threads=2), expected)

    # Trying row-major dense LHS matrix.
    tdptr1 <- initializeCpp(t(DelayedArray(t(d1))))
    expect_equal(tatami.multiply(tdptr1, mat, right=FALSE, num.threads=1), expected)
    expect_equal(tatami.multiply(tdptr1, mat, right=FALSE, num.threads=2), expected)

    # Trying two matrices.
    ptr2 <- initializeCpp(mat)
    expect_equal(tatami.multiply(ptr1, ptr2, right=FALSE, num.threads=1), expected)
    expect_equal(tatami.multiply(ptr1, ptr2, right=FALSE, num.threads=2), expected)
})
