# This tests the ability of the API to properly access character matrices of different types.
# library(testthat); library(beachtest); source("test-character-input.R")

library(beachtest)
genwords <- function(n = 5000) {
    all.choices <- c(rep("", 4), LETTERS) # to get variable length strings.
    a <- do.call(paste0, replicate(5, sample(all.choices, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

#######################################################
# Testing simple matrices:

set.seed(12345)
sFUN <- function(nr=15, nc=10, lambda=5) {
    matrix(genwords(nr*nc), nr, nc)
}

test_that("Simple character matrix input is okay", {
    check_read_all(sFUN, mode="character")
    check_read_all(sFUN, nr=5, nc=30, mode="character")
    check_read_all(sFUN, nr=30, nc=5, mode="character")

    check_read_slice(sFUN, mode="character")
    check_read_slice(sFUN, nr=5, nc=30, mode="character")
    check_read_slice(sFUN, nr=30, nc=5, mode="character")

    check_read_varslice(sFUN, mode="character")
    check_read_varslice(sFUN, nr=5, nc=30, mode="character")
    check_read_varslice(sFUN, nr=30, nc=5, mode="character")

    check_read_const(sFUN, mode="character")
    check_read_const(sFUN, nr=5, nc=30, mode="character")
    check_read_const(sFUN, nr=30, nc=5, mode="character")

    check_read_indexed(sFUN, mode="character")
    check_read_indexed(sFUN, nr=5, nc=30, mode="character")
    check_read_indexed(sFUN, nr=30, nc=5, mode="character")

    check_read_multi(sFUN, mode="character")
    check_read_multi(sFUN, nr=5, nc=30, mode="character")
    check_read_multi(sFUN, nr=30, nc=5, mode="character")

    check_read_type(sFUN, mode="character")
    check_read_class(sFUN(), mode="character", "simple")

    check_read_errors(sFUN, mode="character")
    check_read_all(sFUN, nr=0, nc=0, mode="character")
    check_read_all(sFUN, nr=10, nc=0, mode="character")
    check_read_all(sFUN, nr=0, nc=10, mode="character")
})

#######################################################
# Testing RLE matrices, treated as unknown.

set.seed(23456)
library(DelayedArray)
rFUN <- function(nr=15, nc=10, lambda=1, chunk.ncols=NULL) {
    x <- matrix(sample(LETTERS[1:4], nr*nc, replace=TRUE), nr, nc)
    rle <- Rle(x)
    RleArray(rle, dim(x))
}

test_that("RLE character matrix input is okay", {
    expect_s4_class(rFUN(), "RleMatrix")

    check_read_all(rFUN, mode="character")
    check_read_all(rFUN, nr=5, nc=30, mode="character")
    check_read_all(rFUN, nr=30, nc=5, mode="character")

    check_read_slice(rFUN, mode="character")
    check_read_slice(rFUN, nr=5, nc=30, mode="character")
    check_read_slice(rFUN, nr=30, nc=5, mode="character")

    check_read_varslice(rFUN, mode="character")
    check_read_varslice(rFUN, nr=5, nc=30, mode="character")
    check_read_varslice(rFUN, nr=30, nc=5, mode="character")

    check_read_const(rFUN, mode="character")
    check_read_const(rFUN, nr=5, nc=30, mode="character")
    check_read_const(rFUN, nr=30, nc=5, mode="character")

    check_read_indexed(rFUN, mode="character")
    check_read_indexed(rFUN, nr=5, nc=30, mode="character")
    check_read_indexed(rFUN, nr=30, nc=5, mode="character")

    check_read_multi(rFUN, mode="character")
    check_read_multi(rFUN, nr=5, nc=30, mode="character")
    check_read_multi(rFUN, nr=30, nc=5, mode="character")

    check_read_type(rFUN, mode="character")
    check_read_class(rFUN(), mode="character", "unknown")

    check_read_errors(rFUN, mode="character")
    check_read_all(rFUN, nr=0, nc=0, mode="character")
    check_read_all(rFUN, nr=10, nc=0, mode="character")
    check_read_all(rFUN, nr=0, nc=10, mode="character")
})

test_that("RLE character matrix input is okay with reduced block size", {
    old <- getOption("DelayedArray.block.size")
    options(DelayedArray.block.size=6*8) # seems that characters use 8 byte addresses.

    check_read_all(rFUN, mode="character")
    check_read_all(rFUN, nr=5, nc=30, mode="character")
    check_read_all(rFUN, nr=30, nc=5, mode="character")

    check_read_slice(rFUN, mode="character")
    check_read_slice(rFUN, nr=5, nc=30, mode="character")
    check_read_slice(rFUN, nr=30, nc=5, mode="character")

    check_read_varslice(rFUN, mode="character")
    check_read_varslice(rFUN, nr=5, nc=30, mode="character")
    check_read_varslice(rFUN, nr=30, nc=5, mode="character")

    check_read_const(rFUN, mode="character")
    check_read_const(rFUN, nr=5, nc=30, mode="character")
    check_read_const(rFUN, nr=30, nc=5, mode="character")

    check_read_indexed(rFUN, mode="character")
    check_read_indexed(rFUN, nr=5, nc=30, mode="character")
    check_read_indexed(rFUN, nr=30, nc=5, mode="character")

    check_read_multi(rFUN, mode="character")
    check_read_multi(rFUN, nr=5, nc=30, mode="character")
    check_read_multi(rFUN, nr=30, nc=5, mode="character")

    check_read_type(rFUN, mode="character")
    check_read_class(rFUN(), mode="character", "unknown")

    check_read_errors(rFUN, mode="character")
    check_read_all(rFUN, nr=0, nc=0, mode="character")
    check_read_all(rFUN, nr=10, nc=0, mode="character")
    check_read_all(rFUN, nr=0, nc=10, mode="character")
            
    options(DelayedArray.block.size=old)
})

#######################################################
# Testing HDF5 matrices:

set.seed(34567)
library(HDF5Array)
hFUN <- function(nr=15, nc=10) {
    as(sFUN(nr, nc), "HDF5Array")
}

test_that("HDF5 character matrix input is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_read_all(hFUN, mode="character")
    check_read_all(hFUN, nr=5, nc=30, mode="character")
    check_read_all(hFUN, nr=30, nc=5, mode="character")

    check_read_slice(hFUN, mode="character")
    check_read_slice(hFUN, nr=5, nc=30, mode="character")
    check_read_slice(hFUN, nr=30, nc=5, mode="character")

    check_read_varslice(hFUN, mode="character")
    check_read_varslice(hFUN, nr=5, nc=30, mode="character")
    check_read_varslice(hFUN, nr=30, nc=5, mode="character")

    check_read_const(hFUN, mode="character")
    check_read_const(hFUN, nr=5, nc=30, mode="character")
    check_read_const(hFUN, nr=30, nc=5, mode="character")

    check_read_indexed(hFUN, mode="character")
    check_read_indexed(hFUN, nr=5, nc=30, mode="character")
    check_read_indexed(hFUN, nr=30, nc=5, mode="character")

    check_read_multi(hFUN, mode="character")
    check_read_multi(hFUN, nr=5, nc=30, mode="character")
    check_read_multi(hFUN, nr=30, nc=5, mode="character")

    check_read_type(hFUN, mode="character")
    check_read_class(hFUN(), mode="character", "HDF5")

    check_read_errors(hFUN, mode="character")
    check_read_all(hFUN, nr=0, nc=0, mode="character")
    check_read_all(hFUN, nr=10, nc=0, mode="character")
    check_read_all(hFUN, nr=0, nc=10, mode="character")
})

#######################################################
# Testing delayed operations:

set.seed(981347)
library(DelayedArray)
test_that("Delayed character matrix input is okay", {
    delfuns <- delayed_funs(sFUN, DELAYED_FUN=DelayedArray::tolower)

    for (FUN in delfuns) {
        NR <- 10 + sample(10, 1)
        NC <- 10 + sample(10, 1)
        expect_s4_class(FUN(), "DelayedMatrix")

        check_read_all(FUN, NR, NC, mode="character")

        check_read_slice(FUN, NR, NC, mode="character")

        check_read_varslice(FUN, NR, NC, mode="character")

        check_read_const(FUN, NR, NC, mode="character")

        check_read_indexed(FUN, NR, NC, mode="character")

        check_read_multi(FUN, NR, NC, mode="character")

        check_read_type(FUN, NR, NC, mode="character")
        check_read_class(FUN(), mode="character", "delayed")

        check_read_errors(FUN, NR, NC, mode="character")
        check_read_all(FUN, nr=0, nc=0, mode="character")
        check_read_all(FUN, nr=10, nc=0, mode="character")
        check_read_all(FUN, nr=0, nc=10, mode="character")
    }

    # Proper type check upon coercion!
    expect_identical("logical", .Call("get_type", hFUN()=="A", PACKAGE="beachtest"))
})
