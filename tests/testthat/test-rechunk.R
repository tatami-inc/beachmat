# This checks that the rechunking function is working.

library(HDF5Array)
test_that("rechunking is working", {
    set.seed(1000)

    # Square chunks:
    A <- writeHDF5Array(matrix(runif(5000), nrow=100, ncol=50), chunk=c(10, 10))
    ref <- as.matrix(A)

    byrow <- rechunkByMargins(A, byrow=TRUE)
    expect_identical(ref, as.matrix(byrow))
    byrow <- rechunkByMargins(A, byrow=TRUE, size=23)
    expect_identical(ref, as.matrix(byrow))

    bycol <- rechunkByMargins(A, byrow=FALSE)
    expect_identical(ref, as.matrix(bycol))
    bycol <- rechunkByMargins(A, byrow=FALSE, size=20)
    expect_identical(ref, as.matrix(bycol))

    # Pure column chunks:
    B <- writeHDF5Array(matrix(runif(5000), nrow=100, ncol=50), chunk=c(100, 1))
    ref <- as.matrix(B)

    byrow <- rechunkByMargins(B, byrow=TRUE)
    expect_identical(ref, as.matrix(byrow))
    byrow <- rechunkByMargins(B, byrow=TRUE, size=13)
    expect_identical(ref, as.matrix(byrow))

    bycol <- rechunkByMargins(B, byrow=FALSE)
    expect_identical(ref, as.matrix(bycol))
    bycol <- rechunkByMargins(B, byrow=FALSE, size=19)
    expect_identical(ref, as.matrix(bycol))

    # Pure row chunks:
    C <- writeHDF5Array(matrix(runif(5000), nrow=100, ncol=50), chunk=c(1, 50))
    ref <- as.matrix(C)

    byrow <- rechunkByMargins(C, byrow=TRUE)
    expect_identical(ref, as.matrix(byrow))
    byrow <- rechunkByMargins(C, byrow=TRUE, size=15)
    expect_identical(ref, as.matrix(byrow))

    bycol <- rechunkByMargins(C, byrow=FALSE)
    expect_identical(ref, as.matrix(bycol))
    bycol <- rechunkByMargins(C, byrow=FALSE, size=17)
    expect_identical(ref, as.matrix(bycol))

    # Non-perfect multiple of dimensions and chunk width.
    D <- writeHDF5Array(matrix(runif(5000), nrow=100, ncol=50), chunk=c(21, 11))
    ref <- as.matrix(D)

    byrow <- rechunkByMargins(D, byrow=TRUE)
    expect_identical(ref, as.matrix(byrow))
    byrow <- rechunkByMargins(D, byrow=TRUE, size=11)
    expect_identical(ref, as.matrix(byrow))

    bycol <- rechunkByMargins(D, byrow=FALSE)
    expect_identical(ref, as.matrix(bycol))
    bycol <- rechunkByMargins(D, byrow=FALSE, size=25)
    expect_identical(ref, as.matrix(bycol))

    # Checking for contiguous input.
    E <- writeHDF5Array(matrix(runif(5000), nrow=100, ncol=50), level=0)
    ref <- as.matrix(E)

    byrow <- rechunkByMargins(E, byrow=TRUE)
    expect_identical(ref, as.matrix(byrow))
    byrow <- rechunkByMargins(E, byrow=TRUE, size=20)
    expect_identical(ref, as.matrix(byrow))

    bycol <- rechunkByMargins(E, byrow=FALSE)
    expect_identical(ref, as.matrix(bycol))
    bycol <- rechunkByMargins(E, byrow=FALSE, size=10)
    expect_identical(ref, as.matrix(bycol))

    expect_error(rechunkByMargins(D, outlevel=0), "compression level of 0 implies a contiguous layout")
})
