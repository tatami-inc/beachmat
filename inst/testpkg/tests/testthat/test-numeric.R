# This tests the ability of the API to properly access numeric matrices of different types.
# library(testthat); library(beachtest); source("test-numeric.R")

#######################################################

# Testing simple matrices:

set.seed(12345)
sFUN <- function(nr=15, nc=10) {
    matrix(rnorm(nr*nc), nr, nc)
}

test_that("Simple numeric matrix input is okay", {
    check_numeric_mat(sFUN)
    check_numeric_mat(sFUN, nr=5, nc=30)
    check_numeric_mat(sFUN, nr=30, nc=5)

    check_numeric_slice(sFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))

    # Testing const and non-zero options.   
    check_numeric_const_mat(sFUN)
    check_numeric_const_slice(sFUN, by.row=list(1:5, 6:8))
    
    check_numeric_indexed_mat(sFUN)
    check_numeric_indexed_slice(sFUN, by.row=list(1:5, 6:8))

    check_type(sFUN, expected="double")
})

# Testing dense matrices:

set.seed(13579)
library(Matrix)
dFUN <- function(nr=15, nc=10) {
    Matrix(sFUN(nr, nc), sparse=FALSE, doDiag=FALSE)
}

test_that("Dense numeric matrix input is okay", {
    expect_s4_class(dFUN(), "dgeMatrix")

    check_numeric_mat(dFUN)
    check_numeric_mat(dFUN, nr=5, nc=30)
    check_numeric_mat(dFUN, nr=30, nc=5)

    check_numeric_slice(dFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))

    # Testing const and non-zero options.   
    check_numeric_const_mat(dFUN)
    check_numeric_const_slice(dFUN, by.row=list(1:5, 6:8))
    
    check_numeric_indexed_mat(dFUN)
    check_numeric_indexed_slice(dFUN, by.row=list(1:5, 6:8))

    check_type(dFUN, expected="double")
})

# Testing sparse matrices (dgCMatrix):

set.seed(23456)
csFUN <- function(nr=15, nc=10, d=0.1) {
    rsparsematrix(nrow=nr, ncol=nc, density=d)
}

test_that("Sparse numeric matrix input is okay", {
    expect_s4_class(csFUN(), "dgCMatrix")
    
    check_numeric_mat(csFUN)
    check_numeric_mat(csFUN, nr=5, nc=30)
    check_numeric_mat(csFUN, nr=30, nc=5)
    
    check_numeric_mat(csFUN, d=0.2)
    check_numeric_mat(csFUN, nr=5, nc=30, d=0.2)
    check_numeric_mat(csFUN, nr=30, nc=5, d=0.2)
    
    check_numeric_slice(csFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))

    # Testing const and non-zero options.   
    check_numeric_const_mat(csFUN)
    check_numeric_const_slice(csFUN, by.row=list(1:5, 6:8))
    
    check_numeric_indexed_mat(csFUN)
    check_numeric_indexed_slice(csFUN, by.row=list(1:5, 6:8))
   
    check_type(csFUN, expected="double")
})

# Testing dense symmetric matrices (dspMatrix):

set.seed(45678)
spFUN <- function(nr=10, mode="U") {
    pack(forceSymmetric(sFUN(nr, nr), uplo=mode))
}

test_that("Symmetric numeric matrix input is okay", {
    expect_s4_class(spFUN(), "dspMatrix")

    check_numeric_mat(spFUN)
    check_numeric_mat(spFUN, nr=5)
    check_numeric_mat(spFUN, nr=30)
    
    check_numeric_mat(spFUN, mode="L")
    check_numeric_mat(spFUN, nr=5, mode="L")
    check_numeric_mat(spFUN, nr=30, mode="L")
    
    check_numeric_slice(spFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))
    check_numeric_slice(spFUN, mode="L", by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))
    
    # Testing const and non-zero options.   
    check_numeric_const_mat(spFUN)
    check_numeric_const_mat(spFUN, mode="L")
    check_numeric_const_slice(spFUN, by.row=list(1:5, 6:8))
    check_numeric_const_slice(spFUN, mode="L", by.row=list(1:5, 6:8))
    
    check_numeric_indexed_mat(spFUN)
    check_numeric_indexed_mat(spFUN, mode="L")
    check_numeric_indexed_slice(spFUN, by.row=list(1:5, 6:8))
    check_numeric_indexed_slice(spFUN, mode="L", by.row=list(1:5, 6:8))
 
    check_type(spFUN, expected="double")
    check_type(spFUN, mode="L", expected="double")
})

# Testing Rle matrices:

set.seed(23456)
library(DelayedArray)
rFUN <- function(nr=15, nc=10, density=0.2, chunk.ncols=NULL) {
    x <- as.matrix(csFUN(nr, nc, density))
    rle <- Rle(x)
    if (!is.null(chunk.ncols)) {
        chunksize <- chunk.ncols*nrow(x)
    } else {
        chunksize <- NULL
    }
    RleArray(rle, dim(x), chunksize=chunksize)
}

test_that("RLE numeric matrix input is okay", {
    check_numeric_mat(rFUN)
    check_numeric_mat(rFUN, nr=5, nc=30)
    check_numeric_mat(rFUN, nr=30, nc=5)

    check_numeric_mat(rFUN, density=0.1)
    check_numeric_mat(rFUN, nr=5, nc=30, density=0.1)
    check_numeric_mat(rFUN, nr=30, nc=5, density=0.1)
    
    check_numeric_slice(rFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))
    check_numeric_slice(rFUN, density=0.1, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))

    check_numeric_const_mat(rFUN)
    check_numeric_const_slice(rFUN, by.row=list(1:5, 6:8))
    
    check_numeric_indexed_mat(rFUN)
    check_numeric_indexed_slice(rFUN, by.row=list(1:5, 6:8))

    # Testing chunk settings.
    check_numeric_mat(rFUN, chunk.ncols=3)
    check_numeric_mat(rFUN, nr=5, nc=30, chunk.ncols=5)
    check_numeric_mat(rFUN, nr=30, nc=5, chunk.ncols=2)

    check_numeric_mat(rFUN, density=0.1, chunk.ncols=3)
    check_numeric_mat(rFUN, nr=5, nc=30, density=0.1, chunk.ncols=5)
    check_numeric_mat(rFUN, nr=30, nc=5, density=0.1, chunk.ncols=2)
    
    check_numeric_slice(rFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8), chunk.ncols=2)
    check_numeric_slice(rFUN, density=0.1, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8), chunk.ncols=2)
 
    check_numeric_const_mat(rFUN, chunk.ncols=2)
    check_numeric_const_slice(rFUN, chunk.ncols=2, by.row=list(1:5, 6:8))
    
    check_numeric_indexed_mat(rFUN, chunk.ncols=2)
    check_numeric_indexed_slice(rFUN, chunk.ncols=2, by.row=list(1:5, 6:8))
  
    # Checking type.
    check_type(rFUN, expected="double")
    check_type(rFUN, chunk.ncols=2, expected="double")
})

# Testing HDF5 matrices:

set.seed(34567)
library(HDF5Array)
hFUN <- function(nr=15, nc=10) {
    as(sFUN(nr, nc), "HDF5Array")
}

test_that("HDF5 numeric matrix input is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_numeric_mat(hFUN)
    check_numeric_mat(hFUN, nr=5, nc=30)
    check_numeric_mat(hFUN, nr=30, nc=5)
    
    check_numeric_slice(hFUN, by.row=list(1:5, 6:8), by.col=list(1:5, 6:8))

    # Checking const and non-zero options.
    check_numeric_const_mat(hFUN)
    check_numeric_const_slice(hFUN, by.row=list(1:5, 6:8))
    
    check_numeric_indexed_mat(hFUN)
    check_numeric_indexed_slice(hFUN, by.row=list(1:5, 6:8))

    check_type(hFUN, expected="double")
})

# Testing delayed operations

set.seed(91283)
library(DelayedArray)
test_that("Delayed numeric matrix input is okay", {
    # HDF5-based seed.
    hdf5.funs <- delayed_funs(hFUN)
    for (FUN in hdf5.funs) {
        expect_s4_class(FUN(), "DelayedMatrix")
        check_numeric_mat(FUN)
        check_type(FUN, expected="double")
    }

    # Sparse seed.
    sparse.funs <- delayed_funs(csFUN)
    for (FUN in sparse.funs) {
        expect_s4_class(FUN(), "DelayedMatrix")
        check_numeric_mat(FUN)
        check_type(FUN, expected="double")
    }

    # Simple seed.
    simple.funs <- delayed_funs(sFUN)
    for (FUN in simple.funs) {
        expect_s4_class(FUN(), "DelayedMatrix")
        check_numeric_mat(FUN)
        check_type(FUN, expected="double")
    }
    
    # Proper type check upon coercion!
    expect_identical("logical", .Call(cxx_test_type_check, hFUN() > 0)) 
})

#######################################################

# Testing conversions.

test_that("Numeric matrix input conversions are okay", {
    check_numeric_conversion(sFUN)

    check_numeric_conversion(dFUN)

    check_numeric_conversion(csFUN)

    check_numeric_conversion(spFUN)

    check_numeric_conversion(hFUN)
})

# Testing error generation.

test_that("Numeric matrix input error generation is okay", {
    check_numeric_edge_errors(sFUN)

    check_numeric_edge_errors(dFUN)

    check_numeric_edge_errors(csFUN)

    check_numeric_edge_errors(spFUN)

    check_numeric_edge_errors(hFUN)
})

#######################################################

# Testing simple numeric output:

set.seed(12345)

test_that("Simple numeric matrix output is okay", {
    check_numeric_output_mat(sFUN)
    check_numeric_output_mat(sFUN, nr=5, nc=30)
    
    check_numeric_output_slice(sFUN, by.row=1:12, by.col=3:7)
    check_numeric_output_slice(sFUN, nr=5, nc=30, by.row=2:4, by.col=12:25)
    
    check_numeric_output_indexed(sFUN, N=c(5, 10, 20))
    check_numeric_output_indexed(sFUN, nr=5, nc=30, N=c(5, 10, 20))
})

# Testing HDF5 sparse output:

test_that("Sparse numeric matrix output is okay", {
    check_numeric_output_mat(csFUN)
    check_numeric_output_mat(csFUN, d=0.2)
    check_numeric_output_mat(csFUN, d=0.5)

    check_numeric_output_slice(csFUN, by.row=1:5, by.col=7:9)
    check_numeric_output_slice(csFUN, by.row=1, by.col=2:8, d=0.2)
    check_numeric_output_slice(csFUN, by.row=3:9, by.col=5, d=0.5)

    check_numeric_output_indexed(csFUN, d=0.2, N=c(5, 10, 20))
    check_numeric_output_indexed(csFUN, d=0.5, nr=5, nc=30, N=c(5, 10, 20))
    check_numeric_output_indexed(csFUN, d=0.5, N=c(5, 10, 20))
})

# Testing HDF5 numeric output:

test_that("HDF5 numeric matrix output is okay", {
    check_numeric_output_mat(hFUN)
    check_numeric_output_mat(hFUN, nr=5, nc=30)
    
    check_numeric_output_slice(hFUN, by.row=1:2, by.col=2:10)
    check_numeric_output_slice(hFUN, nr=5, nc=30, by.row=1:2, by.col=2:10)

    check_numeric_output_indexed(hFUN, N=c(5, 10, 20))
    check_numeric_output_indexed(hFUN, nr=5, nc=30, N=c(5, 10, 20))

    check_numeric_order(hFUN)
})

# Testing conversions:

test_that("Numeric matrix output conversions are okay", {
    check_numeric_converted_output(sFUN)
    
    check_numeric_converted_output(hFUN)
})

# Testing mode choices:

test_that("Numeric matrix mode choices are okay", {
    expect_identical(check_output_mode(sFUN, simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(check_output_mode(csFUN, simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(check_output_mode(spFUN, simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(check_output_mode(rFUN, simplify=TRUE, preserve.zero=FALSE), "simple")
    expect_identical(check_output_mode(csFUN, simplify=FALSE, preserve.zero=TRUE), "sparse")
    expect_identical(check_output_mode(csFUN, simplify=FALSE, preserve.zero=FALSE), "HDF5")
    expect_identical(check_output_mode(rFUN, simplify=FALSE, preserve.zero=FALSE), "HDF5")
    expect_identical(check_output_mode(spFUN, simplify=FALSE, preserve.zero=FALSE), "HDF5")
    expect_identical(check_output_mode(hFUN, simplify=FALSE, preserve.zero=FALSE), "HDF5")
})

# Testing for errors:

test_that("Numeric matrix output error generation is okay", {
    check_numeric_edge_output_errors(sFUN)
    
    check_numeric_edge_output_errors(hFUN)
})

#######################################################

