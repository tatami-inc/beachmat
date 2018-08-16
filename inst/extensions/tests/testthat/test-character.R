# Tests whether Aaron's Matrix can be successfully read by beachmat for strings.
# library(testthat); library(morebeach); source("test-character")
#testpkg <- system.file("testpkg", package="beachmat")
#devtools::install(testpkg, quick=TRUE)

library(beachtest); library(morebeach)
generator <- function(nr=15, nc=10) {
    AaronMatrix(matrix(sample(LETTERS, nr*nc, replace=TRUE), nr, nc))
}

test_that("Aaron's Matrix can be read by beachmat", {
    check_read_all(generator, mode="character")
    check_read_all(generator, nr=5, nc=30, mode="character")
    check_read_all(generator, nr=30, nc=5, mode="character")

    check_read_slice(generator, mode="character")
    check_read_slice(generator, nr=5, nc=30, mode="character")
    check_read_slice(generator, nr=30, nc=5, mode="character")

    check_read_varslice(generator, mode="character")
    check_read_varslice(generator, nr=5, nc=30, mode="character")
    check_read_varslice(generator, nr=30, nc=5, mode="character")

    check_read_const(generator, mode="character")
    check_read_const(generator, nr=5, nc=30, mode="character")
    check_read_const(generator, nr=30, nc=5, mode="character")

    check_read_indexed(generator, mode="character")
    check_read_indexed(generator, nr=5, nc=30, mode="character")
    check_read_indexed(generator, nr=30, nc=5, mode="character")

    check_read_multi(generator, mode="character")
    check_read_multi(generator, nr=5, nc=30, mode="character")
    check_read_multi(generator, nr=30, nc=5, mode="character")

    check_read_type(generator, mode="character")
    check_read_class(generator(), mode="character", "external")

    check_read_errors(generator, mode="character")
    check_read_all(generator, nr=0, nc=0, mode="character")
    check_read_all(generator, nr=10, nc=0, mode="character")
    check_read_all(generator, nr=0, nc=10, mode="character")
})
