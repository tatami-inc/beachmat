# This tests the ability of the API to properly output integer matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-integer-output.R")

sFUN <- integer_sFUN
rFUN <- integer_rFUN
hFUN <- integer_hFUN

#######################################################

set.seed(12346)
test_that("Simple integer matrix output is okay", {
    check_write_all(sFUN, mode="integer")
    check_write_all(sFUN, nr=5, nc=30, mode="integer")
    check_write_all(sFUN, nr=30, nc=5, mode="integer")

    check_write_slice(sFUN, mode="integer")
    check_write_slice(sFUN, nr=5, nc=30, mode="integer")
    check_write_slice(sFUN, nr=30, nc=5, mode="integer")

    check_write_varslice(sFUN, mode="integer")
    check_write_varslice(sFUN, nr=5, nc=30, mode="integer")
    check_write_varslice(sFUN, nr=30, nc=5, mode="integer")

    check_write_indexed(sFUN, mode="integer")
    check_write_indexed(sFUN, nr=5, nc=30, mode="integer")
    check_write_indexed(sFUN, nr=30, nc=5, mode="integer")

    check_write_type(sFUN, mode="integer")
    check_write_errors(sFUN, mode="integer")

    check_write_all(sFUN, nr=0, nc=0, mode="integer")
    check_write_all(sFUN, nr=10, nc=0, mode="integer")
    check_write_all(sFUN, nr=0, nc=10, mode="integer")
})

#######################################################

set.seed(34568)
test_that("HDF5 integer matrix output is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_write_all(hFUN, mode="integer")
    check_write_all(hFUN, nr=5, nc=30, mode="integer")
    check_write_all(hFUN, nr=30, nc=5, mode="integer")

    check_write_slice(hFUN, mode="integer")
    check_write_slice(hFUN, nr=5, nc=30, mode="integer")
    check_write_slice(hFUN, nr=30, nc=5, mode="integer")

    check_write_varslice(hFUN, mode="integer")
    check_write_varslice(hFUN, nr=5, nc=30, mode="integer")
    check_write_varslice(hFUN, nr=30, nc=5, mode="integer")

    check_write_indexed(hFUN, mode="integer")
    check_write_indexed(hFUN, nr=5, nc=30, mode="integer")
    check_write_indexed(hFUN, nr=30, nc=5, mode="integer")

    check_write_type(hFUN, mode="integer")
    check_write_errors(hFUN, mode="integer")

    check_write_HDF5(hFUN, mode="integer")

    check_write_all(hFUN, nr=0, nc=0, mode="integer")
    check_write_all(hFUN, nr=10, nc=0, mode="integer")
    check_write_all(hFUN, nr=0, nc=10, mode="integer")
})

#######################################################

test_that("Integer matrix mode choices are okay", {
    check_write_class(sFUN(), "simple", simplify=TRUE)
    check_write_class(sFUN(), "simple", simplify=FALSE)
    check_write_class(sFUN(), "simple", preserve.zeroes=FALSE)

    check_write_class(hFUN(), "HDF5", simplify=TRUE) 
    check_write_class(hFUN(), "HDF5", simplify=FALSE) 

    check_write_class(rFUN(), "simple", simplify=TRUE) 
    check_write_class(rFUN(), "HDF5", simplify=FALSE) 

    check_write_class(hFUN()+1L, "simple", simplify=TRUE) 
    check_write_class(hFUN()+1L, "HDF5", simplify=FALSE) 
})
