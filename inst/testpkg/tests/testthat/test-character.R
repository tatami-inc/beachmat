# This tests the ability of the API to properly access character matrices of different types.
# library(testthat); source("test-character.R")

genwords <- function(n = 5000) {
    all.choices <- c(rep("", 4), LETTERS) # to get variable length strings.
    a <- do.call(paste0, replicate(5, sample(all.choices, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

#######################################################

# Testing simple matrices:

set.seed(12345)
sFUN <- function(nr=15, nc=10) {
    matrix(genwords(nr*nc), nr, nc)
}

test_that("Simple character matrix input is okay", {
    beachtest:::check_character_mat(sFUN)
    beachtest:::check_character_mat(sFUN, nr=5, nc=30)
    beachtest:::check_character_mat(sFUN, nr=30, nc=5)
    
    beachtest:::check_character_slice(sFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))

    # Testing const options.
    beachtest:::check_character_const_mat(sFUN)
    beachtest:::check_character_const_slice(sFUN, by.row=list(1:5, 6:8))

    beachtest:::check_type(sFUN, expected="character")
})

# Testing RLE matrices:

set.seed(23456)
library(DelayedArray)
rFUN <- function(nr=15, nc=10, chunk.ncols=NULL) {
    x <- matrix(sample(LETTERS[1:4], nr*nc, replace=TRUE), nr, nc)
    rle <- Rle(x)
    if (!is.null(chunk.ncols)) {
        chunksize <- chunk.ncols*nrow(x)
    } else {
        chunksize <- NULL
    }
    RleArray(rle, dim(x), chunksize=chunksize)
}

test_that("RLE character matrix input is okay", {
    beachtest:::check_character_mat(rFUN)
    beachtest:::check_character_mat(rFUN, nr=5, nc=30)
    beachtest:::check_character_mat(rFUN, nr=30, nc=5)
    
    beachtest:::check_character_slice(rFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))

    beachtest:::check_character_const_mat(rFUN)
    beachtest:::check_character_const_slice(rFUN, by.row=list(1:5, 6:8))

    # Testing chunk settings.
    beachtest:::check_character_mat(rFUN, chunk.ncol=3)
    beachtest:::check_character_mat(rFUN, nr=5, nc=30, chunk.ncol=5)
    beachtest:::check_character_mat(rFUN, nr=30, nc=5, chunk.ncol=2)
    
    beachtest:::check_character_slice(rFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8), chunk.ncol=2)

    beachtest:::check_character_const_mat(rFUN, chunk.ncols=2)
    beachtest:::check_character_const_slice(rFUN, chunk.ncols=2, by.row=list(1:5, 6:8))

    # Checking type. 
    beachtest:::check_type(rFUN, expected="character")
    beachtest:::check_type(rFUN, chunk.ncol=2, expected="character")
})

# Testing HDF5 matrices:

set.seed(34567)
library(HDF5Array)
hFUN <- function(nr=15, nc=10) {
    as(sFUN(nr, nc), "HDF5Array")
}

test_that("HDF5 character matrix input is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")
    beachtest:::check_character_mat(hFUN)
    beachtest:::check_character_mat(hFUN, nr=5, nc=30)
    beachtest:::check_character_mat(hFUN, nr=30, nc=5)
    
    beachtest:::check_character_slice(hFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))
    
    # Testing const options.
    beachtest:::check_character_const_mat(hFUN)
    beachtest:::check_character_const_slice(hFUN, by.row=list(1:5, 6:8))

    beachtest:::check_type(hFUN, expected="character")
})

# Testing delayed operations:

set.seed(981347)
library(DelayedArray)
test_that("Delayed character matrix input is okay", {
    # HDF5-based seed.
    hdf5.funs <- beachtest:::delayed_funs(hFUN)
    for (FUN in hdf5.funs) {
        expect_s4_class(FUN(), "DelayedMatrix")
        beachtest:::check_character_mat(FUN)
        beachtest:::check_type(FUN, expected="character")
    }

    # Simple seed.
    simple.funs <- beachtest:::delayed_funs(sFUN)
    for (FUN in simple.funs) {
        expect_s4_class(FUN(), "DelayedMatrix")
        beachtest:::check_character_mat(FUN)
        beachtest:::check_type(FUN, expected="character")
    }

    # Trigger realization.
    alt_hFUN <- function(..., transpose=FALSE) {
        out <- DelayedArray::tolower(hFUN(...))
        if (transpose) { 
            out <- DelayedArray::t(out)
        }
        return(out)
    }
    expect_s4_class(alt_hFUN(), "DelayedMatrix")
    beachtest:::check_character_mat(alt_hFUN)
    beachtest:::check_type(alt_hFUN, expected="character")
 
    expect_s4_class(alt_hFUN(transpose=TRUE), "DelayedMatrix")
    beachtest:::check_character_mat(alt_hFUN, transpose=TRUE) # checking that transposition WITH delayed ops wipes the transformer.
    beachtest:::check_type(alt_hFUN, transpose=TRUE, expected="character")

    comb_hFUN <- function(...) {
        DelayedArray::cbind(hFUN(...), hFUN(...))
    }
    expect_s4_class(comb_hFUN(), "DelayedMatrix")
    beachtest:::check_character_mat(comb_hFUN) # checking that odd seed types are properly realized.
    beachtest:::check_type(comb_hFUN, expected="character")
     
    # Proper type check!
    expect_identical("logical", .Call(beachtest:::cxx_test_type_check, hFUN()=="A")) 
})

#######################################################

# Testing errors.

test_that("Character matrix input error generation is okay", {
    beachtest:::check_character_edge_errors(sFUN)
    beachtest:::check_character_edge_errors(hFUN)
})

#######################################################

# Testing simple character output:

set.seed(12345)

test_that("Simple character matrix output is okay", {
    beachtest:::check_character_output_mat(sFUN)
    beachtest:::check_character_output_mat(sFUN, nr=5, nc=30)

    beachtest:::check_character_output_slice(sFUN, by.row=2:11, by.col=4:8)
    beachtest:::check_character_output_slice(sFUN, nr=5, nc=30, by.row=2:3, by.col=4:28)
})

# Testing HDF5 character output:

test_that("HDF5 character matrix output is okay", {
    beachtest:::check_character_output_mat(hFUN)
    beachtest:::check_character_output_mat(hFUN, nr=5, nc=30)

    beachtest:::check_character_output_slice(hFUN, by.row=5:15, by.col=8:10)
    beachtest:::check_character_output_slice(hFUN, nr=5, nc=30, by.row=2:5, by.col=8:30)

    beachtest:::check_character_order(hFUN)
})

test_that("Numeric matrix mode choices are okay", {
    expect_identical(beachtest:::check_output_mode(sFUN, simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(beachtest:::check_output_mode(rFUN, simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(beachtest:::check_output_mode(rFUN, simplify=FALSE, preserve.zero=FALSE), "HDF5")
    expect_identical(beachtest:::check_output_mode(hFUN, simplify=FALSE, preserve.zero=FALSE), "HDF5")
})

# Testing for errors:

test_that("Character matrix output error generation is okay", {
    beachtest:::check_character_edge_output_errors(sFUN)

    beachtest:::check_character_edge_output_errors(hFUN)
})

#######################################################

