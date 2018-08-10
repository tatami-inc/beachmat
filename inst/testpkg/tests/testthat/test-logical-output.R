# This tests the ability of the API to properly output logical matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-logical-output.R")

sFUN <- logical_sFUN
dFUN <- logical_dFUN
csFUN <- logical_csFUN
tsFUN <- logical_tsFUN
hFUN <- logical_hFUN

#######################################################

set.seed(12346)
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

    check_write_all(sFUN, nr=0, nc=0, mode="logical")
    check_write_all(sFUN, nr=10, nc=0, mode="logical")
    check_write_all(sFUN, nr=0, nc=10, mode="logical")
})

#######################################################

set.seed(23457)
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

    check_write_type(csFUN, mode="logical")
    check_write_errors(csFUN, mode="logical")

    check_write_all(csFUN, nr=0, nc=0, mode="logical")
    check_write_all(csFUN, nr=10, nc=0, mode="logical")
    check_write_all(csFUN, nr=0, nc=10, mode="logical")
})

#######################################################

set.seed(34568)
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

    check_write_all(hFUN, nr=0, nc=0, mode="logical")
    check_write_all(hFUN, nr=10, nc=0, mode="logical")
    check_write_all(hFUN, nr=0, nc=10, mode="logical")

    check_write_HDF5(hFUN, mode="logical")
})

#######################################################

test_that("Logical matrix mode choices are okay", {
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

    check_write_class(!hFUN(), "simple", simplify=TRUE) 
    check_write_class(!hFUN(), "HDF5", simplify=FALSE) 
})
