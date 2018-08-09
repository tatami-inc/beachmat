# This tests the ability of the API to properly output logical matrices of different types.
# library(testthat); library(beachtest); source("test-logical-output.R")

library(beachtest)

#######################################################
# Testing simple matrices:

set.seed(12345)
sFUN <- function(nr=15, nc=10) {
    matrix(rbinom(nr*nc, 1, 0.5)==0, nr, nc)
}

test_that("Simple logical matrix output is okay", {
    check_write_all(sFUN, mode="logical")
    check_write_all(sFUN, nr=5, nc=30, mode="logical")
    check_write_all(sFUN, nr=30, nc=5, mode="logical")

    check_write_slice(sFUN, mode="logical")
    check_write_slice(sFUN, nr=5, nc=30, mode="logical")
    check_write_slice(sFUN, nr=30, nc=5, mode="logical")

    check_write_varslice(sFUN, mode="logical")
    check_write_varslice(sFUN, nr=5, nc=30, mode="logical")
    check_write_varslice(sFUN, nr=30, nc=5, mode="logical")

    check_write_indexed(sFUN, mode="logical")
    check_write_indexed(sFUN, nr=5, nc=30, mode="logical")
    check_write_indexed(sFUN, nr=30, nc=5, mode="logical")

    check_write_type(sFUN, mode="logical")
    check_write_errors(sFUN, mode="logical")
})

#######################################################
# Testing sparse matrices.

set.seed(23456)
library(Matrix)
csFUN <- function(nr=15, nc=10, d=0.1) {
    rsparsematrix(nrow=nr, ncol=nc, density=d) != 0
}

test_that("sparse logical matrix output is okay", {
    expect_s4_class(csFUN(), "lgCMatrix")

    check_write_all(csFUN, mode="logical")
    check_write_all(csFUN, nr=5, nc=30, mode="logical")
    check_write_all(csFUN, nr=30, nc=5, mode="logical")

    check_write_slice(csFUN, mode="logical")
    check_write_slice(csFUN, nr=5, nc=30, mode="logical")
    check_write_slice(csFUN, nr=30, nc=5, mode="logical")

    check_write_varslice(csFUN, mode="logical")
    check_write_varslice(csFUN, nr=5, nc=30, mode="logical")
    check_write_varslice(csFUN, nr=30, nc=5, mode="logical")

    check_write_indexed(csFUN, mode="logical")
    check_write_indexed(csFUN, nr=5, nc=30, mode="logical")
    check_write_indexed(csFUN, nr=30, nc=5, mode="logical")

    check_write_errors(csFUN, mode="logical")

    # Not checking type conversions, as integers don't get sparse support and logical->logical isn't consistent between C++ and R.
})

#######################################################
# Testing RLE matrices, treated as unknown.

set.seed(23456)
library(DelayedArray)
rFUN <- function(nr=15, nc=10) {
    x <- sFUN(nr, nc)
    rle <- Rle(x)
    RleArray(rle, dim(x))
}

test_that("RLE logical matrix output is okay", {
    expect_s4_class(rFUN(), "RleMatrix")

    check_write_all(rFUN, mode="logical", out.class="matrix")
    check_write_all(rFUN, nr=5, nc=30, mode="logical", out.class="matrix")
    check_write_all(rFUN, nr=30, nc=5, mode="logical", out.class="matrix")

    check_write_slice(rFUN, mode="logical", out.class="matrix")
    check_write_slice(rFUN, nr=5, nc=30, mode="logical", out.class="matrix")
    check_write_slice(rFUN, nr=30, nc=5, mode="logical", out.class="matrix")

    check_write_varslice(rFUN, mode="logical", out.class="matrix")
    check_write_varslice(rFUN, nr=5, nc=30, mode="logical", out.class="matrix")
    check_write_varslice(rFUN, nr=30, nc=5, mode="logical", out.class="matrix")

    check_write_indexed(rFUN, mode="logical", out.class="matrix")
    check_write_indexed(rFUN, nr=5, nc=30, mode="logical", out.class="matrix")
    check_write_indexed(rFUN, nr=30, nc=5, mode="logical", out.class="matrix")

    check_write_type(rFUN, mode="logical", out.class="matrix")
    check_write_errors(rFUN, mode="logical", out.class="matrix")
})

#######################################################
# Testing HDF5 matrices:

set.seed(34567)
library(HDF5Array)
hFUN <- function(nr=15, nc=10) {
    as(sFUN(nr, nc), "HDF5Array")
}

test_that("HDF5 logical matrix output is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_write_all(hFUN, mode="logical")
    check_write_all(hFUN, nr=5, nc=30, mode="logical")
    check_write_all(hFUN, nr=30, nc=5, mode="logical")

    check_write_slice(hFUN, mode="logical")
    check_write_slice(hFUN, nr=5, nc=30, mode="logical")
    check_write_slice(hFUN, nr=30, nc=5, mode="logical")

    check_write_varslice(hFUN, mode="logical")
    check_write_varslice(hFUN, nr=5, nc=30, mode="logical")
    check_write_varslice(hFUN, nr=30, nc=5, mode="logical")

    check_write_indexed(hFUN, mode="logical")
    check_write_indexed(hFUN, nr=5, nc=30, mode="logical")
    check_write_indexed(hFUN, nr=30, nc=5, mode="logical")

    check_write_type(hFUN, mode="logical")
    check_write_errors(hFUN, mode="logical")

    check_write_HDF5(hFUN, mode="logical")
})
