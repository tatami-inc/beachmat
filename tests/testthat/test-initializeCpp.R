# This checks the initialization procedure.
# library(testthat); library(beachmat); source("test-initializeCpp.R")

am_i_ok <- function(ref, ptr, exact=TRUE) {
    expect_identical(dim(ref), beachmat:::tatami_dim(ptr))
    test <- if (exact) expect_identical else expect_equal
    for (i in seq_len(ncol(ref))) {
        expected <- ref[,i]
        if (!is.double(expected)) { 
            expected <- as.double(expected) 
        }
        test(expected, beachmat:::tatami_column(ptr, i))
    }
}

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

    de <- Matrix(dd)
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

test_that("initialization works correctly with scalar arithmetic", {
    z0 <- DelayedArray(y)

    z <- z0 + 1
    ptr <- initializeCpp(z)
    am_i_ok(y + 1, ptr)

    z <- z0 * 10
    ptr <- initializeCpp(z)
    am_i_ok(y * 10, ptr)

    z <- z0 - 1
    ptr <- initializeCpp(z)
    am_i_ok(y - 1, ptr)

    z <- 1 - z0
    ptr <- initializeCpp(z)
    am_i_ok(1 - y, ptr)

    z <- z0 / 10
    ptr <- initializeCpp(z)
    am_i_ok(y / 10, ptr)

    z <- 10 / z0
    ptr <- initializeCpp(z)
    am_i_ok(10 / y, ptr)
})

test_that("initialization works correctly with vector arithmetic", {
    z0 <- DelayedArray(y)

    {
        vr <- runif(nrow(y))

        z <- z0 + vr
        ptr <- initializeCpp(z)
        am_i_ok(y + vr, ptr)

        z <- z0 * vr
        ptr <- initializeCpp(z)
        am_i_ok(y * vr, ptr)

        z <- z0 - vr
        ptr <- initializeCpp(z)
        am_i_ok(y - vr, ptr)

        z <- vr - z0
        ptr <- initializeCpp(z)
        am_i_ok(vr - y, ptr)

        z <- z0 / vr
        ptr <- initializeCpp(z)
        am_i_ok(y / vr, ptr)

        z <- vr / z0
        ptr <- initializeCpp(z)
        am_i_ok(vr / y, ptr)
    }

    {
        vc <- runif(ncol(y))

        z <- sweep(z0, 2, vc, "+")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) + vc), ptr)

        z <- sweep(z0, 2, vc, "*")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) * vc), ptr)

        z <- sweep(z0, 2, vc, "-")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) - vc), ptr)

        z <- sweep(z0, 2, vc, "/")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) / vc), ptr)
    }
})

test_that("initialization works correctly with DelayedArray log-transforms", {
    z0 <- DelayedArray(y)

    # Adding a value to make sure it's not log-zero.
    z <- log2(z0 + 2)
    ptr <- initializeCpp(z)
    am_i_ok(log2(y + 2), ptr, exact=FALSE)

    z <- log10(z0 + 2)
    ptr <- initializeCpp(z)
    am_i_ok(log10(y + 2), ptr, exact=FALSE)

    z <- log(z0 + 2)
    ptr <- initializeCpp(z)
    am_i_ok(log(y + 2), ptr, exact=FALSE)

    z <- log(z0 + 2, 2.5)
    ptr <- initializeCpp(z)
    am_i_ok(log(y + 2, 2.5), ptr, exact=FALSE)
})

test_that("initialization works correctly with DelayedArray rounding", {
    ref <- (DelayedArray(x) + 0.001) * 12.3 # avoid problems at 0.5.

    x0 <- DelayedArray(ref)
    z <- round(x0)
    ptr <- initializeCpp(z)
    am_i_ok(round(ref), ptr)

    z <- round(x0, 2)
    expect_error(initializeCpp(z), "digits = 0")
})

test_that("initialization works correctly with other DelayedArray unary operations", {
    z0 <- DelayedArray(y)

    z <- +z0 
    ptr <- initializeCpp(z)
    am_i_ok(y, ptr)

    z <- -z0 
    ptr <- initializeCpp(z)
    am_i_ok(-y, ptr)

    z <- log1p(z0)
    ptr <- initializeCpp(z)
    am_i_ok(log1p(y), ptr, exact=FALSE)

    z <- sqrt(z0)
    ptr <- initializeCpp(z)
    am_i_ok(sqrt(y), ptr, exact=FALSE)

    z <- exp(z0)
    ptr <- initializeCpp(z)
    am_i_ok(exp(y), ptr, exact=FALSE)

    x0 <- DelayedArray(x)
    z <- abs(x0)
    ptr <- initializeCpp(z)
    am_i_ok(abs(x), ptr)
})
