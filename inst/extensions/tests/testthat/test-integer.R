# Tests whether Aaron's Matrix can be successfully read by beachmat.
# library(testthat); library(morebeach); source("test-integer")
#testpkg <- system.file("testpkg", package="beachmat")
#devtools::install(testpkg, quick=TRUE)

library(beachtest); library(morebeach)
generator <- function(nr=15, nc=10) {
    AaronMatrix(matrix(rpois(nr*nc, lambda=5), nr, nc))
}

test_that("Aaron's Matrix can be read by beachmat", {
    check_read_all(generator, mode="integer")
    check_read_all(generator, nr=5, nc=30, mode="integer")
    check_read_all(generator, nr=30, nc=5, mode="integer")

    check_read_slice(generator, mode="integer")
    check_read_slice(generator, nr=5, nc=30, mode="integer")
    check_read_slice(generator, nr=30, nc=5, mode="integer")

    check_read_varslice(generator, mode="integer")
    check_read_varslice(generator, nr=5, nc=30, mode="integer")
    check_read_varslice(generator, nr=30, nc=5, mode="integer")

    check_read_const(generator, mode="integer")
    check_read_const(generator, nr=5, nc=30, mode="integer")
    check_read_const(generator, nr=30, nc=5, mode="integer")

    check_read_indexed(generator, mode="integer")
    check_read_indexed(generator, nr=5, nc=30, mode="integer")
    check_read_indexed(generator, nr=30, nc=5, mode="integer")

    check_read_multi(generator, mode="integer")
    check_read_multi(generator, nr=5, nc=30, mode="integer")
    check_read_multi(generator, nr=30, nc=5, mode="integer")

    check_read_type(generator, mode="integer")
#    check_read_class(generator(), mode="integer", "external")

#   check_read_errors(generator, mode="integer")
    check_read_all(generator, nr=0, nc=0, mode="integer")
    check_read_all(generator, nr=10, nc=0, mode="integer")
    check_read_all(generator, nr=0, nc=10, mode="integer")
})
