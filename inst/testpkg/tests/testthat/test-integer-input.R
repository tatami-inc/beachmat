# This tests the ability of the API to properly access integer matrices of different types.
# library(testthat); library(beachtest); source("test-integer-input.R")

library(beachtest)

#######################################################
# Testing simple matrices:

set.seed(12345)
sFUN <- function(nr=15, nc=10, lambda=5) {
    matrix(rpois(nr*nc, lambda=lambda), nr, nc)
}

test_that("Simple integer matrix input is okay", {
    check_read_all(sFUN, mode="integer")
    check_read_all(sFUN, nr=5, nc=30, mode="integer")
    check_read_all(sFUN, nr=30, nc=5, mode="integer")

    check_read_slice(sFUN, mode="integer")
    check_read_slice(sFUN, nr=5, nc=30, mode="integer")
    check_read_slice(sFUN, nr=30, nc=5, mode="integer")

    check_read_varslice(sFUN, mode="integer")
    check_read_varslice(sFUN, nr=5, nc=30, mode="integer")
    check_read_varslice(sFUN, nr=30, nc=5, mode="integer")

    check_read_const(sFUN, mode="integer")
    check_read_const(sFUN, nr=5, nc=30, mode="integer")
    check_read_const(sFUN, nr=30, nc=5, mode="integer")

    check_read_indexed(sFUN, mode="integer")
    check_read_indexed(sFUN, nr=5, nc=30, mode="integer")
    check_read_indexed(sFUN, nr=30, nc=5, mode="integer")

    check_read_multi(sFUN, mode="integer")
    check_read_multi(sFUN, nr=5, nc=30, mode="integer")
    check_read_multi(sFUN, nr=30, nc=5, mode="integer")

    check_read_type(sFUN, mode="integer")
    check_read_errors(sFUN, mode="integer")

    check_read_all(sFUN, nr=0, nc=0, mode="integer")
    check_read_all(sFUN, nr=10, nc=0, mode="integer")
    check_read_all(sFUN, nr=0, nc=10, mode="integer")
})

#######################################################
# Testing RLE matrices, treated as unknown.

set.seed(23456)
library(DelayedArray)
rFUN <- function(nr=15, nc=10, lambda=1, chunk.ncols=NULL) {
    x <- sFUN(nr, nc, lambda=lambda)
    rle <- Rle(x)
    RleArray(rle, dim(x))
}

test_that("RLE integer matrix input is okay", {
    expect_s4_class(rFUN(), "RleMatrix")

    check_read_all(rFUN, mode="integer")
    check_read_all(rFUN, nr=5, nc=30, mode="integer")
    check_read_all(rFUN, nr=30, nc=5, mode="integer")

    check_read_slice(rFUN, mode="integer")
    check_read_slice(rFUN, nr=5, nc=30, mode="integer")
    check_read_slice(rFUN, nr=30, nc=5, mode="integer")

    check_read_varslice(rFUN, mode="integer")
    check_read_varslice(rFUN, nr=5, nc=30, mode="integer")
    check_read_varslice(rFUN, nr=30, nc=5, mode="integer")

    check_read_const(rFUN, mode="integer")
    check_read_const(rFUN, nr=5, nc=30, mode="integer")
    check_read_const(rFUN, nr=30, nc=5, mode="integer")

    check_read_indexed(rFUN, mode="integer")
    check_read_indexed(rFUN, nr=5, nc=30, mode="integer")
    check_read_indexed(rFUN, nr=30, nc=5, mode="integer")

    check_read_multi(rFUN, mode="integer")
    check_read_multi(rFUN, nr=5, nc=30, mode="integer")
    check_read_multi(rFUN, nr=30, nc=5, mode="integer")

    check_read_type(rFUN, mode="integer")
    check_read_errors(rFUN, mode="integer")

    check_read_all(rFUN, nr=0, nc=0, mode="integer")
    check_read_all(rFUN, nr=10, nc=0, mode="integer")
    check_read_all(rFUN, nr=0, nc=10, mode="integer")
})

test_that("RLE integer matrix input is okay with reduced block size", {
    old <- getOption("DelayedArray.block.size")
    options(DelayedArray.block.size=6*4)

    check_read_all(rFUN, mode="integer")
    check_read_all(rFUN, nr=5, nc=30, mode="integer")
    check_read_all(rFUN, nr=30, nc=5, mode="integer")

    check_read_slice(rFUN, mode="integer")
    check_read_slice(rFUN, nr=5, nc=30, mode="integer")
    check_read_slice(rFUN, nr=30, nc=5, mode="integer")

    check_read_varslice(rFUN, mode="integer")
    check_read_varslice(rFUN, nr=5, nc=30, mode="integer")
    check_read_varslice(rFUN, nr=30, nc=5, mode="integer")

    check_read_const(rFUN, mode="integer")
    check_read_const(rFUN, nr=5, nc=30, mode="integer")
    check_read_const(rFUN, nr=30, nc=5, mode="integer")

    check_read_indexed(rFUN, mode="integer")
    check_read_indexed(rFUN, nr=5, nc=30, mode="integer")
    check_read_indexed(rFUN, nr=30, nc=5, mode="integer")

    check_read_multi(rFUN, mode="integer")
    check_read_multi(rFUN, nr=5, nc=30, mode="integer")
    check_read_multi(rFUN, nr=30, nc=5, mode="integer")

    check_read_type(rFUN, mode="integer")
    check_read_errors(rFUN, mode="integer")
             
    check_read_all(rFUN, nr=0, nc=0, mode="integer")
    check_read_all(rFUN, nr=10, nc=0, mode="integer")
    check_read_all(rFUN, nr=0, nc=10, mode="integer")

    options(DelayedArray.block.size=old)
})

#######################################################
# Testing HDF5 matrices:

set.seed(34567)
library(HDF5Array)
hFUN <- function(nr=15, nc=10) {
    as(sFUN(nr, nc), "HDF5Array")
}

test_that("HDF5 integer matrix input is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_read_all(hFUN, mode="integer")
    check_read_all(hFUN, nr=5, nc=30, mode="integer")
    check_read_all(hFUN, nr=30, nc=5, mode="integer")

    check_read_slice(hFUN, mode="integer")
    check_read_slice(hFUN, nr=5, nc=30, mode="integer")
    check_read_slice(hFUN, nr=30, nc=5, mode="integer")

    check_read_varslice(hFUN, mode="integer")
    check_read_varslice(hFUN, nr=5, nc=30, mode="integer")
    check_read_varslice(hFUN, nr=30, nc=5, mode="integer")

    check_read_const(hFUN, mode="integer")
    check_read_const(hFUN, nr=5, nc=30, mode="integer")
    check_read_const(hFUN, nr=30, nc=5, mode="integer")

    check_read_indexed(hFUN, mode="integer")
    check_read_indexed(hFUN, nr=5, nc=30, mode="integer")
    check_read_indexed(hFUN, nr=30, nc=5, mode="integer")

    check_read_multi(hFUN, mode="integer")
    check_read_multi(hFUN, nr=5, nc=30, mode="integer")
    check_read_multi(hFUN, nr=30, nc=5, mode="integer")

    check_read_type(hFUN, mode="integer")
    check_read_errors(hFUN, mode="integer")

    check_read_all(hFUN, nr=0, nc=0, mode="integer")
    check_read_all(hFUN, nr=10, nc=0, mode="integer")
    check_read_all(hFUN, nr=0, nc=10, mode="integer")
})

#######################################################
# Testing delayed operations:

set.seed(981347)
library(DelayedArray)
test_that("Delayed integer matrix input is okay", {
    delfuns <- delayed_funs(sFUN, DELAYED_FUN=function(x) { x + sample(nrow(x)) })

    for (FUN in delfuns) {
        NR <- 10 + sample(10, 1)
        NC <- 10 + sample(10, 1)
        expect_s4_class(FUN(), "DelayedMatrix")

        check_read_all(FUN, NR, NC, mode="integer")

        check_read_slice(FUN, NR, NC, mode="integer")

        check_read_varslice(FUN, NR, NC, mode="integer")

        check_read_const(FUN, NR, NC, mode="integer")

        check_read_indexed(FUN, NR, NC, mode="integer")

        check_read_multi(FUN, NR, NC, mode="integer")

        check_read_type(FUN, NR, NC, mode="integer")

        check_read_errors(FUN, NR, NC, mode="integer")

        # Edge case checks.
        check_read_all(FUN, nr=0, nc=0, mode="integer")
        check_read_all(FUN, nr=10, nc=0, mode="integer")
        check_read_all(FUN, nr=0, nc=10, mode="integer")
    }

    # Proper type check upon coercion!
    expect_identical("double", .Call("get_type", hFUN()+1, PACKAGE="beachtest"))
})
