# This checks the initialization procedure.
# library(testthat); library(beachmat); source("setup.R"); source("test-initializeCpp-isometric_binary.R")

library(DelayedArray)
set.seed(1000)
x1 <- Matrix::rsparsematrix(1000, 100, 0.1)
y1 <- round(abs(x1)*10)
x2 <- Matrix::rsparsematrix(1000, 100, 0.1)
y2 <- round(abs(x2)*10)

test_that("initialization works correctly with binary arithmetic", {
    z01 <- DelayedArray(y1)
    z02 <- DelayedArray(y2)

    z <- z01 + z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 + y2, ptr)

    z <- z01 * z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 * y2, ptr)

    z <- z01 - z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 - y2, ptr)

    z <- z01 / z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 / y2, ptr)

    z <- z01 ^ z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 ^ y2, ptr, exact=(.Platform$OS.type != "windows")) # I don't make the rules.

    z <- z01 %% z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 %% y2, ptr)

    z <- z01 %/% z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 %/% y2, ptr)
})

test_that("initialization works correctly with binary comparisons", {
    z01 <- DelayedArray(y1)
    z02 <- DelayedArray(y2)

    z <- z01 > z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 > y2, ptr)

    z <- z01 >= z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 >= y2, ptr)

    z <- z01 < z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 < y2, ptr)

    z <- z01 <= z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 <= y2, ptr)

    z <- z01 == z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 == y2, ptr)

    z <- z01 != z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 != y2, ptr)
})

test_that("initialization works correctly with binary boolean operations", {
    z01 <- DelayedArray(y1)
    z02 <- DelayedArray(y2)

    z <- z01 | z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 | y2, ptr)

    z <- z01 & z02
    ptr <- initializeCpp(z)
    am_i_ok(y1 & y2, ptr)
})
