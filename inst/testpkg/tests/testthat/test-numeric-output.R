# This tests the ability of the API to properly output numeric matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-numeric-output.R")

sFUN <- numeric_sFUN
dFUN <- numeric_dFUN
csFUN <- numeric_csFUN
tsFUN <- numeric_tsFUN
hFUN <- numeric_hFUN

#######################################################

set.seed(12346)
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

set.seed(23457)
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

    check_write_type(csFUN, mode="numeric")
    check_write_errors(csFUN, mode="numeric")

    check_write_all(csFUN, nr=0, nc=0, mode="numeric")
    check_write_all(csFUN, nr=10, nc=0, mode="numeric")
    check_write_all(csFUN, nr=0, nc=10, mode="numeric")
})

#######################################################

set.seed(34568)
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
    check_write_class(sFUN(), "simple", simplify=TRUE)
    check_write_class(sFUN(), "simple", simplify=FALSE)
    check_write_class(sFUN(), "simple", preserve.zeroes=FALSE)

    check_write_class(dFUN(), "simple", simplify=TRUE)
    check_write_class(dFUN(), "simple", simplify=FALSE)
    check_write_class(dFUN(), "simple", preserve.zeroes=FALSE)

    check_write_class(csFUN(), "simple", simplify=TRUE, preserve.zeroes=FALSE) 
    check_write_class(csFUN(), "HDF5", simplify=FALSE, preserve.zeroes=FALSE) 
    check_write_class(csFUN(), "sparse", simplify=FALSE, preserve.zeroes=TRUE) 
    check_write_class(csFUN(), "sparse", simplify=FALSE, preserve.zeroes=TRUE) 

    check_write_class(hFUN(), "HDF5", simplify=TRUE) 
    check_write_class(hFUN(), "HDF5", simplify=FALSE) 
    
    check_write_class(tsFUN(), "simple", simplify=TRUE) 
    check_write_class(tsFUN(), "HDF5", simplify=FALSE) 

    check_write_class(hFUN()+1, "simple", simplify=TRUE) 
    check_write_class(hFUN()+1, "HDF5", simplify=FALSE) 
})
