# This checks the initialization procedure.
# library(testthat); library(beachmat); source("setup.R"); source("test-initializeCpp-other.R")

set.seed(1000)
x <- Matrix::rsparsematrix(1000, 100, 0.1)
y <- round(abs(x)*10)

test_that("initialization works correctly with dense matrices", {
    dd <- as.matrix(y)
    {
        ptr <- initializeCpp(dd)
        am_i_ok(dd, ptr)

        dd2 <- dd
        storage.mode(dd2) <- "integer"
        ptr <- initializeCpp(dd2)
        am_i_ok(dd2, ptr)

        dd2 <- dd
        storage.mode(dd2) <- "logical"
        ptr <- initializeCpp(dd2)
        am_i_ok(dd2, ptr)
    }

    de <- Matrix::Matrix(dd)
    {
        ptr <- initializeCpp(de)
        am_i_ok(de, ptr)

        de2 <- de > 0
        ptr <- initializeCpp(de2)
        am_i_ok(de2, ptr)
    }
})

test_that("initialization works correctly with sparse matrices", {
    {
        ptr <- initializeCpp(y)
        am_i_ok(y, ptr)

        z <- new("dgRMatrix", x=y@x, j=y@i, p=y@p, Dim=rev(y@Dim))
        ptr <- initializeCpp(z)
        am_i_ok(z, ptr)
    }

    {
        y2 <- y != 0
        ptr <- initializeCpp(y2)
        am_i_ok(y2, ptr)

        z2 <- new("lgRMatrix", x=y2@x, j=y2@i, p=y2@p, Dim=rev(y2@Dim))
        ptr <- initializeCpp(z2)
        am_i_ok(z2, ptr)
    }
})

test_that("initialization works correctly with SVT sparse matrices", {
    library(SparseArray)

    {
        z <- as(y, "SVT_SparseMatrix")
        ptr <- initializeCpp(z)
        am_i_ok(y, ptr)
    }

    {
        y2 <- y != 0
        z <- as(y2, "SVT_SparseMatrix")
        ptr <- initializeCpp(y2)
        am_i_ok(y2, ptr)
    }

    {
        y2 <- as.matrix(y)
        storage.mode(y2) <- "integer"
        z <- as(y2, "SVT_SparseMatrix")
        ptr <- initializeCpp(y2)
        am_i_ok(y2, ptr)
    }
})

library(DelayedArray)
test_that("initialization works correctly with DelayedArray", {
    z <- DelayedArray(y)
    ptr <- initializeCpp(z)
    am_i_ok(y, ptr)
})

test_that("initialization works correctly with DelayedArray transposition", {
    z0 <- DelayedArray(y)
    z <- t(z0)
    ptr <- initializeCpp(z)
    am_i_ok(t(y), ptr)
})

test_that("initialization works correctly with DelayedArray subsetting", {
    z0 <- DelayedArray(y)

    rkeep <- sample(nrow(y), 100)
    z <- z0[rkeep,]
    ptr <- initializeCpp(z)
    am_i_ok(y[rkeep,], ptr)

    ckeep <- sample(ncol(y), 10)
    z <- z0[,ckeep]
    ptr <- initializeCpp(z)
    am_i_ok(y[,ckeep], ptr)

    z <- z0[rkeep,ckeep]
    ptr <- initializeCpp(z)
    am_i_ok(y[rkeep,ckeep], ptr)

    rkeep <- 100:200
    ckeep <- 5:30
    z <- z0[rkeep,ckeep]
    ptr <- initializeCpp(z)
    am_i_ok(y[rkeep,ckeep], ptr)
})

test_that("initialization works correctly with DelayedArray combining", {
    z0 <- DelayedArray(y)

    x2 <- Matrix::rsparsematrix(99, 100, 0.1)
    y2 <- round(abs(x)*10)
    z <- rbind(z0, DelayedArray(y2))
    ptr <- initializeCpp(z)
    am_i_ok(rbind(y, y2), ptr)

    x2 <- Matrix::rsparsematrix(1000, 50, 0.1)
    y2 <- round(abs(x)*10)
    z <- cbind(z0, DelayedArray(y2))
    ptr <- initializeCpp(z)
    am_i_ok(cbind(y, y2), ptr)
})

test_that("initialization works correctly with an unknown DelayedArray", {
    library(HDF5Array)
    mat <- matrix(rnorm(50), ncol=5)
    mat2 <- as(mat, "HDF5Array")

    ptr <- initializeCpp(mat2)
    am_i_ok(mat, ptr)

    # works in the sparse case.
    mat <- DelayedArray(Matrix::rsparsematrix(100, 50, 0.1))
    mat2 <- round(mat, digits=2)

    expect_warning(ptr <- initializeCpp(mat2), "falling back")
    am_i_ok(mat2, ptr)
})

test_that("initialization no-ops correctly with its own output", {
    dd <- as.matrix(y)
    ptr <- initializeCpp(dd)
    ptr2 <- initializeCpp(ptr)
    am_i_ok(dd, ptr2)
})
