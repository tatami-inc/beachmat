# This tests the ability of the API to properly output character matrices of different types.
# library(testthat); library(beachtest); source("test-character-output.R")

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

test_that("Simple character matrix output is okay", {
    check_write_all(sFUN, mode="character")
    check_write_all(sFUN, nr=5, nc=30, mode="character")
    check_write_all(sFUN, nr=30, nc=5, mode="character")

    check_write_slice(sFUN, mode="character")
    check_write_slice(sFUN, nr=5, nc=30, mode="character")
    check_write_slice(sFUN, nr=30, nc=5, mode="character")

    check_write_varslice(sFUN, mode="character")
    check_write_varslice(sFUN, nr=5, nc=30, mode="character")
    check_write_varslice(sFUN, nr=30, nc=5, mode="character")

    check_write_indexed(sFUN, mode="character")
    check_write_indexed(sFUN, nr=5, nc=30, mode="character")
    check_write_indexed(sFUN, nr=30, nc=5, mode="character")

    check_write_type(sFUN, mode="character")
    check_write_errors(sFUN, mode="character")

    check_write_all(sFUN, nr=0, nc=0, mode="character")
    check_write_all(sFUN, nr=10, nc=0, mode="character")
    check_write_all(sFUN, nr=0, nc=10, mode="character")
})

#######################################################
# Testing RLE matrices, treated as unknown.

set.seed(23456)
library(DelayedArray)
rFUN <- function(nr=15, nc=10, lambda=1) {
    x <- sFUN(nr, nc, lambda=lambda)
    rle <- Rle(x)
    RleArray(rle, dim(x))
}

test_that("RLE character matrix output is okay", {
    expect_s4_class(rFUN(), "RleMatrix")

    check_write_all(rFUN, mode="character", out.class="matrix")
    check_write_all(rFUN, nr=5, nc=30, mode="character", out.class="matrix")
    check_write_all(rFUN, nr=30, nc=5, mode="character", out.class="matrix")

    check_write_slice(rFUN, mode="character", out.class="matrix")
    check_write_slice(rFUN, nr=5, nc=30, mode="character", out.class="matrix")
    check_write_slice(rFUN, nr=30, nc=5, mode="character", out.class="matrix")

    check_write_varslice(rFUN, mode="character", out.class="matrix")
    check_write_varslice(rFUN, nr=5, nc=30, mode="character", out.class="matrix")
    check_write_varslice(rFUN, nr=30, nc=5, mode="character", out.class="matrix")

    check_write_indexed(rFUN, mode="character", out.class="matrix")
    check_write_indexed(rFUN, nr=5, nc=30, mode="character", out.class="matrix")
    check_write_indexed(rFUN, nr=30, nc=5, mode="character", out.class="matrix")

    check_write_type(rFUN, mode="character", out.class="matrix")
    check_write_errors(rFUN, mode="character", out.class="matrix")

    check_write_all(rFUN, nr=0, nc=0, mode="character", out.class="matrix")
    check_write_all(rFUN, nr=10, nc=0, mode="character", out.class="matrix")
    check_write_all(rFUN, nr=0, nc=10, mode="character", out.class="matrix")
})

#######################################################
# Testing HDF5 matrices:

set.seed(34567)
library(HDF5Array)
hFUN <- function(nr=15, nc=10) {
    as(sFUN(nr, nc), "HDF5Array")
}

test_that("HDF5 character matrix output is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_write_all(hFUN, mode="character")
    check_write_all(hFUN, nr=5, nc=30, mode="character")
    check_write_all(hFUN, nr=30, nc=5, mode="character")

    check_write_slice(hFUN, mode="character")
    check_write_slice(hFUN, nr=5, nc=30, mode="character")
    check_write_slice(hFUN, nr=30, nc=5, mode="character")

    check_write_varslice(hFUN, mode="character")
    check_write_varslice(hFUN, nr=5, nc=30, mode="character")
    check_write_varslice(hFUN, nr=30, nc=5, mode="character")

    check_write_indexed(hFUN, mode="character")
    check_write_indexed(hFUN, nr=5, nc=30, mode="character")
    check_write_indexed(hFUN, nr=30, nc=5, mode="character")

    check_write_type(hFUN, mode="character")
    check_write_errors(hFUN, mode="character")

    check_write_all(hFUN, nr=0, nc=0, mode="character")
    check_write_all(hFUN, nr=10, nc=0, mode="character")
    check_write_all(hFUN, nr=0, nc=10, mode="character")

    check_write_HDF5(hFUN, mode="character")
})

#######################################################

test_that("Character matrix mode choices are okay", {
    check_write_class(sFUN(), "simple", simplify=TRUE)
    check_write_class(sFUN(), "simple", simplify=FALSE)
    check_write_class(sFUN(), "simple", preserve.zeroes=FALSE)

    check_write_class(rFUN(), "simple", simplify=TRUE) 
    check_write_class(rFUN(), "HDF5", simplify=FALSE) 

    check_write_class(hFUN(), "HDF5", simplify=TRUE) 
    check_write_class(hFUN(), "HDF5", simplify=FALSE) 

    check_write_class(DelayedArray::tolower(hFUN()), "simple", simplify=TRUE) 
    check_write_class(DelayedArray::tolower(hFUN()), "HDF5", simplify=FALSE) 
})
