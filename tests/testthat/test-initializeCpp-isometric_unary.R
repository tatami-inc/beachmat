# This checks the initialization procedure.
# library(testthat); library(beachmat); source("setup.R"); source("test-initializeCpp-isometric_unary.R")

set.seed(1000)
x <- Matrix::rsparsematrix(1000, 100, 0.1)
y <- round(abs(x)*10)

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

test_that("initialization works correctly with scalar comparisons", {
    z0 <- DelayedArray(y)

    z <- z0 > 1
    ptr <- initializeCpp(z)
    am_i_ok(y > 1, ptr)

    z <- z0 >= 10
    ptr <- initializeCpp(z)
    am_i_ok(y >= 10, ptr)

    z <- z0 < 3
    ptr <- initializeCpp(z)
    am_i_ok(y < 3, ptr)

    z <- z0 <= 5
    ptr <- initializeCpp(z)
    am_i_ok(y <= 5, ptr)

    z <- z0 == 0
    ptr <- initializeCpp(z)
    am_i_ok(y == 0, ptr)

    z <- z0 != 0
    ptr <- initializeCpp(z)
    am_i_ok(y != 0, ptr)

    # Behaves correctly on the right as well.
    z <- 5 > z0
    ptr <- initializeCpp(z)
    am_i_ok(y < 5, ptr)

    z <- 20 >= z0
    ptr <- initializeCpp(z)
    am_i_ok(y <= 20, ptr)

    z <- 3 < z0
    ptr <- initializeCpp(z)
    am_i_ok(y > 3, ptr)

    z <- 7 <= z0
    ptr <- initializeCpp(z)
    am_i_ok(y >= 7, ptr)

    z <- 0 == z0
    ptr <- initializeCpp(z)
    am_i_ok(y == 0, ptr)

    z <- 0 != z0
    ptr <- initializeCpp(z)
    am_i_ok(y != 0, ptr)
})

test_that("initialization works correctly with vector comparisons", {
    z0 <- DelayedArray(y)

    # By row with the non-array on the left.
    {
        vr <- runif(nrow(y))

        z <- vr > z0 
        ptr <- initializeCpp(z)
        am_i_ok(y < vr, ptr)

        z <- vr >= z0 
        ptr <- initializeCpp(z)
        am_i_ok(y <= vr, ptr)

        z <- vr < z0 
        ptr <- initializeCpp(z)
        am_i_ok(y > vr, ptr)

        z <- vr <= z0
        ptr <- initializeCpp(z)
        am_i_ok(y >= vr, ptr)

        vd <- diag(y)
        z <- vd != z0 
        ptr <- initializeCpp(z)
        am_i_ok(y != vd, ptr)

        z <- vd == z0
        ptr <- initializeCpp(z)
        am_i_ok(vd == y, ptr)
    }

    # By column with the non-array on the right.
    {
        vc <- runif(ncol(y))

        z <- sweep(z0, 2, vc, ">")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) > vc), ptr)

        z <- sweep(z0, 2, vc, ">=")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) >= vc), ptr)

        z <- sweep(z0, 2, vc, "<")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) < vc), ptr)

        z <- sweep(z0, 2, vc, "<=")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) <= vc), ptr)

        vd <- diag(y)
        z <- sweep(z0, 2, vd, "!=")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) != vd), ptr)

        z <- sweep(z0, 2, vd, "==")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) == vd), ptr)
    }
})

test_that("initialization works correctly with scalar boolean operations", {
    z0 <- DelayedArray(y)

    z <- z0 | TRUE
    ptr <- initializeCpp(z)
    am_i_ok(y | TRUE, ptr)

    z <- FALSE & z0
    ptr <- initializeCpp(z)
    am_i_ok(y & FALSE, ptr)

    z <- !z0 
    ptr <- initializeCpp(z)
    am_i_ok(!y, ptr)
})

test_that("initialization works correctly with vector boolean operations", {
    z0 <- DelayedArray(y)

    # By row with the non-array on the left.
    {
        vr <- rbinom(nrow(y), 1, 0.5)

        z <- vr | z0
        ptr <- initializeCpp(z)
        am_i_ok(y | vr, ptr)

        z <- vr & z0
        ptr <- initializeCpp(z)
        am_i_ok(y & vr, ptr)
    }

    # By column with the non-array on the right.
    {
        vc <- rbinom(ncol(y), 1, 0.5)

        z <- sweep(z0, 2, vc, "&")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) & vc), ptr)

        z <- sweep(z0, 2, vc, "|")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) | vc), ptr)
    }
})

test_that("initialization works correctly with DelayedArray exponential and log-transforms", {
    z0 <- DelayedArray(y)

    z <- exp(z0)
    ptr <- initializeCpp(z)
    am_i_ok(exp(y), ptr, exact=FALSE)

    z <- expm1(z0)
    ptr <- initializeCpp(z)
    am_i_ok(expm1(y), ptr, exact=FALSE)

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

    z <- log1p(z0)
    ptr <- initializeCpp(z)
    am_i_ok(log1p(y), ptr, exact=FALSE)
})

test_that("initialization works correctly with DelayedArray rounding", {
    ref <- (DelayedArray(x) + 0.001) * 12.3 # avoid problems at 0.5.

    x0 <- DelayedArray(ref)
    z <- round(x0)
    ptr <- initializeCpp(z)
    am_i_ok(round(ref), ptr)

    z <- round(x0, 2)
    expect_warning(ptr <- initializeCpp(z), "digits = 0")
    am_i_ok(round(ref, 2), ptr)
})

test_that("initialization works correctly with DelayedArray nearest integer operations", {
    x0 <- DelayedArray(x)

    z <- ceiling(x0)
    ptr <- initializeCpp(z)
    am_i_ok(ceiling(x), ptr)

    z <- floor(x0)
    ptr <- initializeCpp(z)
    am_i_ok(floor(x), ptr)

    z <- trunc(x0)
    ptr <- initializeCpp(z)
    am_i_ok(trunc(x), ptr)
})

test_that("initialization works correctly with DelayedArray trigonometric operations", {
    set.seed(1000)
    a <- Matrix::rsparsematrix(1000, 100, 0.1, rand.x = function(n) signif(runif(n, min = -1, max = 1), 2))
    a0 <- DelayedArray(a)

    z <- asin(a0)
    ptr <- initializeCpp(z)
    am_i_ok(asin(a), ptr, exact=FALSE)

    z <- acos(a0)
    ptr <- initializeCpp(z)
    am_i_ok(acos(a), ptr, exact=FALSE)
})

test_that("initialization works correctly with DelayedArray hyperbolic operations", {
    x0 <- DelayedArray(x)

    set.seed(1000)
    a <- Matrix::rsparsematrix(1000, 100, 0.1, rand.x = function(n) signif(runif(n, min = -1, max = 1), 2))
    a0 <- DelayedArray(a)

    z <- asinh(x0)
    ptr <- initializeCpp(z)
    am_i_ok(asinh(x), ptr, exact=FALSE)

    # Adding a value to make sure domain >= 1.
    z <- acosh(a0 + 2)
    ptr <- initializeCpp(z)
    am_i_ok(acosh(a + 2), ptr, exact=FALSE)
})

test_that("initialization works correctly with other DelayedArray unary operations", {
    z0 <- DelayedArray(y)

    z <- +z0 
    ptr <- initializeCpp(z)
    am_i_ok(y, ptr)

    z <- -z0 
    ptr <- initializeCpp(z)
    am_i_ok(-y, ptr)

    z <- sqrt(z0)
    ptr <- initializeCpp(z)
    am_i_ok(sqrt(y), ptr, exact=FALSE)

    x0 <- DelayedArray(x)
    z <- abs(x0)
    ptr <- initializeCpp(z)
    am_i_ok(abs(x), ptr)
})
