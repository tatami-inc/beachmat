# This checks the initialization procedure.
# library(testthat); library(beachmat); source("setup.R"); source("test-initializeCpp-isometric_unary.R")

library(DelayedArray)
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

    z <- z0 ^ 2
    ptr <- initializeCpp(z)
    am_i_ok(y ^ 2, ptr)

    z <- 2 ^ z0
    ptr <- initializeCpp(z)
    am_i_ok(2 ^ y, ptr)

    z <- z0 %% 3
    ptr <- initializeCpp(z)
    am_i_ok(y %% 3, ptr)

    z <- 3 %% z0
    ptr <- initializeCpp(z)
    am_i_ok(3 %% y, ptr)

    z <- z0 %/% 7
    ptr <- initializeCpp(z)
    am_i_ok(y %/% 7, ptr)

    z <- 11 %/% z0
    ptr <- initializeCpp(z)
    am_i_ok(11 %/% y, ptr)
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

        z <- z0 ^ vr
        ptr <- initializeCpp(z)
        am_i_ok(y ^ vr, ptr, exact=FALSE)

        z <- vr ^ z0
        ptr <- initializeCpp(z)
        am_i_ok(vr ^ y, ptr, exact=FALSE)

        z <- z0 %% vr
        ptr <- initializeCpp(z)
        am_i_ok(y %% vr, ptr)

        z <- vr %% z0
        ptr <- initializeCpp(z)
        am_i_ok(vr %% y, ptr)

        z <- z0 %/% vr
        ptr <- initializeCpp(z)
        am_i_ok(y %/% vr, ptr)

        z <- vr %/% z0
        ptr <- initializeCpp(z)
        am_i_ok(vr %/% y, ptr)
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

        z <- sweep(z0, 2, vc, "%%")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) %% vc), ptr)

        z <- sweep(z0, 2, vc, "%/%")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) %/% vc), ptr)

        z <- sweep(z0, 2, vc, "^")
        ptr <- initializeCpp(z)
        am_i_ok(t(t(y) ^ vc), ptr, exact=FALSE)
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

test_that("initialization works correctly with type transformations", {
    x <- Matrix::rsparsematrix(1000, 100, 0.1)
    y <- as.matrix(round(abs(x)*10))

    for (from in c("logical", "integer", "double")) {
        copy <- y
        storage.mode(copy) <- from

        # We don't have native support for coercion to raw... whatever.
        for (to in c("logical", "integer", "double")) {
            z <- DelayedArray(copy)
            type(z) <- to
            expect_warning(ptr <- initializeCpp(z), NA)
            am_i_ok(z, ptr)
        }
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
    expect_warning(ptr <- initializeCpp(z), "digits = 0")
    am_i_ok(round(ref, 2), ptr)
})

test_that("initialization works correctly with other DelayedArray unary operations", {
    y0 <- DelayedArray(y)

    z <- +y0 
    ptr <- initializeCpp(z)
    am_i_ok(y, ptr)

    z <- -y0 
    ptr <- initializeCpp(z)
    am_i_ok(-y, ptr)

    x0 <- DelayedArray(x)
    z <- abs(x0)
    ptr <- initializeCpp(z)
    am_i_ok(abs(x), ptr)

    z <- sign(x0)
    ptr <- initializeCpp(z)
    am_i_ok(sign(x), ptr)

    z <- sqrt(y0)
    ptr <- initializeCpp(z)
    am_i_ok(sqrt(y), ptr, exact=FALSE)

    z <- floor(x0)
    ptr <- initializeCpp(z)
    am_i_ok(floor(x), ptr)

    z <- ceiling(x0)
    ptr <- initializeCpp(z)
    am_i_ok(ceiling(x), ptr)

    z <- trunc(x0)
    ptr <- initializeCpp(z)
    am_i_ok(trunc(x), ptr)

    z <- exp(y0)
    ptr <- initializeCpp(z)
    am_i_ok(exp(y), ptr, exact=FALSE)

    z <- expm1(x0)
    ptr <- initializeCpp(z)
    am_i_ok(expm1(x), ptr, exact=FALSE)

    z <- log1p(y0)
    ptr <- initializeCpp(z)
    am_i_ok(log1p(y), ptr, exact=FALSE)
})

test_that("initialization works correctly with DelayedArray trig operations", {
    # Don't want to deal with different handling of domain/pole errors.
    scaled <- x / (max(abs(x)) + 0.01)
    scaled0 <- DelayedArray(scaled)

    z <- cos(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(cos(scaled), ptr, exact=FALSE)

    z <- sin(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(sin(scaled), ptr, exact=FALSE)

    z <- tan(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(tan(scaled), ptr, exact=FALSE)

    z <- cospi(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(cospi(scaled), ptr, exact=FALSE)

    z <- sinpi(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(sinpi(scaled), ptr, exact=FALSE)

    z <- tanpi(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(tanpi(scaled), ptr, exact=FALSE)

    z <- acos(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(acos(scaled), ptr, exact=FALSE)

    z <- asin(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(asin(scaled), ptr, exact=FALSE)

    z <- atan(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(atan(scaled), ptr, exact=FALSE)

    z <- cosh(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(cosh(scaled), ptr, exact=FALSE)

    z <- sinh(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(sinh(scaled), ptr, exact=FALSE)

    z <- tanh(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(tanh(scaled), ptr, exact=FALSE)

    # Again, no pole errors.
    positive <- DelayedArray(as.matrix(y + 1))

    z <- acosh(positive)
    ptr <- initializeCpp(z)
    am_i_ok(acosh(y + 1), ptr, exact=FALSE)

    z <- asinh(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(asinh(scaled), ptr, exact=FALSE)

    z <- atanh(scaled0)
    ptr <- initializeCpp(z)
    am_i_ok(atanh(scaled), ptr, exact=FALSE)
})

test_that("initialization works correctly with DelayedArray gamma operations", {
    # Don't want to deal with pole errors.
    positive <- DelayedArray(as.matrix(y + 1))

    z <- gamma(positive)
    ptr <- initializeCpp(z)
    am_i_ok(gamma(y + 1), ptr, exact=FALSE)

    z <- lgamma(positive)
    ptr <- initializeCpp(z)
    am_i_ok(lgamma(y + 1), ptr, exact=FALSE)
})
