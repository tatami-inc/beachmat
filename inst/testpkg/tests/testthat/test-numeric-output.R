# This tests the ability of the API to properly output numeric matrices of different types.
# library(testthat); library(beachtest); source("test-numeric-output.R")

library(beachtest)

#######################################################
# Testing simple matrices:

set.seed(12345)
sFUN <- function(nr=15, nc=10) {
    matrix(rnorm(nr*nc), nr, nc)
}

test_that("Simple numeric matrix output is okay", {
    check_write_all(sFUN, mode="numeric")
    check_write_all(sFUN, nr=5, nc=30, mode="numeric")
    check_write_all(sFUN, nr=30, nc=5, mode="numeric")

    check_write_slice(sFUN, mode="numeric")
    check_write_slice(sFUN, nr=5, nc=30, mode="numeric")
    check_write_slice(sFUN, nr=30, nc=5, mode="numeric")

    check_write_varslice(sFUN, mode="numeric")
    check_write_varslice(sFUN, nr=5, nc=30, mode="numeric")
    check_write_varslice(sFUN, nr=30, nc=5, mode="numeric")

    check_write_indexed(sFUN, mode="numeric")
    check_write_indexed(sFUN, nr=5, nc=30, mode="numeric")
    check_write_indexed(sFUN, nr=30, nc=5, mode="numeric")

    check_write_type(sFUN, mode="numeric")
    check_write_errors(sFUN, mode="numeric")

    check_write_all(sFUN, nr=0, nc=0, mode="numeric")
    check_write_all(sFUN, nr=10, nc=0, mode="numeric")
    check_write_all(sFUN, nr=0, nc=10, mode="numeric")
})

#######################################################
# Testing sparse matrices.

set.seed(23456)
library(Matrix)
csFUN <- function(nr=15, nc=10, d=0.1) {
    rsparsematrix(nrow=nr, ncol=nc, density=d)
}

test_that("sparse numeric matrix output is okay", {
    expect_s4_class(csFUN(), "dgCMatrix")

    check_write_all(csFUN, mode="numeric")
    check_write_all(csFUN, nr=5, nc=30, mode="numeric")
    check_write_all(csFUN, nr=30, nc=5, mode="numeric")

    check_write_slice(csFUN, mode="numeric")
    check_write_slice(csFUN, nr=5, nc=30, mode="numeric")
    check_write_slice(csFUN, nr=30, nc=5, mode="numeric")

    check_write_varslice(csFUN, mode="numeric")
    check_write_varslice(csFUN, nr=5, nc=30, mode="numeric")
    check_write_varslice(csFUN, nr=30, nc=5, mode="numeric")

    check_write_indexed(csFUN, mode="numeric")
    check_write_indexed(csFUN, nr=5, nc=30, mode="numeric")
    check_write_indexed(csFUN, nr=30, nc=5, mode="numeric")

    # Not checking type conversions, as integers don't get sparse support and numeric->logical isn't consistent between C++ and R.

    check_write_errors(csFUN, mode="numeric")

    check_write_all(csFUN, nr=0, nc=0, mode="numeric")
    check_write_all(csFUN, nr=10, nc=0, mode="numeric")
    check_write_all(csFUN, nr=0, nc=10, mode="numeric")
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

test_that("RLE numeric matrix output is okay", {
    expect_s4_class(rFUN(), "RleMatrix")

    check_write_all(rFUN, mode="numeric", out.class="matrix")
    check_write_all(rFUN, nr=5, nc=30, mode="numeric", out.class="matrix")
    check_write_all(rFUN, nr=30, nc=5, mode="numeric", out.class="matrix")

    check_write_slice(rFUN, mode="numeric", out.class="matrix")
    check_write_slice(rFUN, nr=5, nc=30, mode="numeric", out.class="matrix")
    check_write_slice(rFUN, nr=30, nc=5, mode="numeric", out.class="matrix")

    check_write_varslice(rFUN, mode="numeric", out.class="matrix")
    check_write_varslice(rFUN, nr=5, nc=30, mode="numeric", out.class="matrix")
    check_write_varslice(rFUN, nr=30, nc=5, mode="numeric", out.class="matrix")

    check_write_indexed(rFUN, mode="numeric", out.class="matrix")
    check_write_indexed(rFUN, nr=5, nc=30, mode="numeric", out.class="matrix")
    check_write_indexed(rFUN, nr=30, nc=5, mode="numeric", out.class="matrix")

    check_write_type(rFUN, mode="numeric", out.class="matrix")
    check_write_errors(rFUN, mode="numeric", out.class="matrix")

    check_write_all(rFUN, nr=0, nc=0, mode="numeric", out.class="matrix")
    check_write_all(rFUN, nr=10, nc=0, mode="numeric", out.class="matrix")
    check_write_all(rFUN, nr=0, nc=10, mode="numeric", out.class="matrix")
})

#######################################################
# Testing HDF5 matrices:

set.seed(34567)
library(HDF5Array)
hFUN <- function(nr=15, nc=10) {
    as(sFUN(nr, nc), "HDF5Array")
}

test_that("HDF5 numeric matrix output is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_write_all(hFUN, mode="numeric")
    check_write_all(hFUN, nr=5, nc=30, mode="numeric")
    check_write_all(hFUN, nr=30, nc=5, mode="numeric")

    check_write_slice(hFUN, mode="numeric")
    check_write_slice(hFUN, nr=5, nc=30, mode="numeric")
    check_write_slice(hFUN, nr=30, nc=5, mode="numeric")

    check_write_varslice(hFUN, mode="numeric")
    check_write_varslice(hFUN, nr=5, nc=30, mode="numeric")
    check_write_varslice(hFUN, nr=30, nc=5, mode="numeric")

    check_write_indexed(hFUN, mode="numeric")
    check_write_indexed(hFUN, nr=5, nc=30, mode="numeric")
    check_write_indexed(hFUN, nr=30, nc=5, mode="numeric")

    check_write_type(hFUN, mode="numeric")
    check_write_errors(hFUN, mode="numeric")

    check_write_all(hFUN, nr=0, nc=0, mode="numeric")
    check_write_all(hFUN, nr=10, nc=0, mode="numeric")
    check_write_all(hFUN, nr=0, nc=10, mode="numeric")

    check_write_HDF5(hFUN, mode="numeric")
})

#######################################################

test_that("Numeric matrix mode choices are okay", {
    check_write_mode(sFUN(), "simple", simplify=TRUE)
    check_write_mode(sFUN(), "simple", simplify=FALSE)
    check_write_mode(sFUN(), "simple", preserve.zeroes=FALSE)

    check_write_mode(csFUN(), "simple", simplify=TRUE, preserve.zeroes=FALSE) 
    check_write_mode(csFUN(), "HDF5", simplify=FALSE, preserve.zeroes=FALSE) 
    check_write_mode(csFUN(), "sparse", simplify=FALSE, preserve.zeroes=TRUE) 
    check_write_mode(csFUN(), "sparse", simplify=FALSE, preserve.zeroes=TRUE) 

    check_write_mode(rFUN(), "simple", simplify=TRUE) 
    check_write_mode(rFUN(), "HDF5", simplify=FALSE) 

    check_write_mode(hFUN(), "HDF5", simplify=TRUE) 
    check_write_mode(hFUN(), "HDF5", simplify=FALSE) 

    check_write_mode(hFUN()+1, "simple", simplify=TRUE) 
    check_write_mode(hFUN()+1, "HDF5", simplify=FALSE) 
})
