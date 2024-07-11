# This tests the blockApply capabilities for a number of different matrices.
# library(testthat); library(beachmat); source("test-apply.R")

library(DelayedArray)
library(BiocParallel)

rs <- function(x, mult=1) Matrix::rowSums(x) * mult

cs <- function(x, mult=1) Matrix::colSums(x) * mult

test_that("apply works with ordinary matrices", {
    x <- matrix(runif(10000), ncol=10)

    # Only one matrix emitted.
    out <- colBlockApply(x, cs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], cs(x))

    out <- rowBlockApply(x, rs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], rs(x))

    # Works with additional arguments.
    out <- colBlockApply(x, cs, mult=2)
    expect_identical(out[[1]], cs(x, 2))

    out <- rowBlockApply(x, rs, mult=2)
    expect_identical(out[[1]], rs(x, 2))
})

test_that("apply on ordinary matrices respects grid construction", {
    x <- matrix(runif(10000), ncol=10)

    # Ignores block size limits.
    setAutoBlockSize(100)

    out <- colBlockApply(x, cs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], cs(x))

    out <- rowBlockApply(x, rs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], rs(x))

    setAutoBlockSize()

    # Works with parallelization.
    BPPARAM <- SnowParam(2)

    out <- colBlockApply(x, cs, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_identical(unlist(out), cs(x))

    out <- rowBlockApply(x, rs, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_identical(unlist(out), rs(x))

    # Parallelization and grid size play nice.
    setAutoBlockSize(nrow(x) * 8)

    out <- colBlockApply(x, cs, grid=TRUE, BPPARAM=BPPARAM)
    expect_identical(length(out), ncol(x))
    expect_identical(unlist(out), cs(x))

    setAutoBlockSize(ncol(x) * 8 * 10)

    out <- rowBlockApply(x, rs, grid=TRUE, BPPARAM=BPPARAM)
    expect_identical(length(out), as.integer(nrow(x) / 10L))
    expect_identical(unlist(out), rs(x))

    setAutoBlockSize()
})

test_that("apply works with sparse matrices", {
    x <- Matrix::rsparsematrix(100, 50, density=0.1)

    # Only one matrix emitted.
    out <- colBlockApply(x, cs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], cs(x))

    out <- rowBlockApply(x, rs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], rs(x))

    # Works with additional arguments.
    out <- colBlockApply(x, cs, mult=2)
    expect_identical(out[[1]], cs(x, 2))

    out <- rowBlockApply(x, rs, mult=2)
    expect_identical(out[[1]], rs(x, 2))
})

test_that("apply on sparse matrices respects grid construction", {
    x <- Matrix::rsparsematrix(100, 50, density=0.1)

    # Ignores block size limits.
    setAutoBlockSize(100)

    out <- colBlockApply(x, cs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], cs(x))

    out <- rowBlockApply(x, rs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], rs(x))

    setAutoBlockSize()

    # Works with parallelization.
    BPPARAM <- SnowParam(2)

    out <- colBlockApply(x, cs, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_identical(unlist(out), cs(x))

    out <- rowBlockApply(x, rs, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_identical(unlist(out), rs(x))

    # Parallelization and grid size play nice.
    setAutoBlockSize(nrow(x) * 8)

    out <- colBlockApply(x, cs, grid=TRUE, BPPARAM=BPPARAM)
    expect_identical(length(out), ncol(x))
    expect_identical(unlist(out), cs(x))

    setAutoBlockSize(ncol(x) * 8 * 10)

    out <- rowBlockApply(x, rs, grid=TRUE, BPPARAM=BPPARAM)
    expect_identical(length(out), as.integer(nrow(x) / 10L))
    expect_identical(unlist(out), rs(x))

    setAutoBlockSize()
})

test_that("apply works with pristine DelayedMatrices", {
    x <- DelayedArray(matrix(runif(10000), ncol=10))

    # Only one matrix emitted.
    out <- colBlockApply(x, cs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], cs(x))

    out <- rowBlockApply(x, rs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], rs(x))

    # Works with additional arguments.
    out <- colBlockApply(x, cs, mult=2)
    expect_identical(out[[1]], cs(x, 2))

    out <- rowBlockApply(x, rs, mult=2)
    expect_identical(out[[1]], rs(x, 2))

    # Works with other backends.
    x <- RleArray(Rle(sample(3, 10000, replace=TRUE)), dim=c(1000, 10))
    out <- colBlockApply(x, cs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], cs(x))

    out <- rowBlockApply(x, rs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], rs(x))
})

test_that("apply on pristine DelayedMatrices respects grid construction", {
    x <- RleArray(Rle(sample(3, 10000, replace=TRUE)), dim=c(1000, 10))

    # Works with parallelization.
    BPPARAM <- SnowParam(2)

    out <- colBlockApply(x, cs, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_identical(unlist(out), cs(x))

    out <- rowBlockApply(x, rs, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_identical(unlist(out), rs(x))

    # Parallelization and grid size play nice.
    setAutoBlockSize(nrow(x) * 4)

    out <- colBlockApply(x, cs, BPPARAM=BPPARAM)
    expect_identical(length(out), ncol(x))
    expect_identical(unlist(out), cs(x))

    setAutoBlockSize(ncol(x) * 4 * 10)

    out <- rowBlockApply(x, rs, BPPARAM=BPPARAM)
    expect_identical(length(out), as.integer(nrow(x) / 10L))
    expect_identical(unlist(out), rs(x))

    setAutoBlockSize()
})

test_that("apply works with non-pristine DelayedMatrices", {
    x <- DelayedArray(matrix(runif(10000), ncol=10)) * 2

    # Only one matrix emitted.
    out <- colBlockApply(x, cs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], cs(x))

    out <- rowBlockApply(x, rs)
    expect_identical(length(out), 1L)
    expect_identical(out[[1]], rs(x))

    # Works with additional arguments.
    out <- colBlockApply(x, cs, mult=2)
    expect_identical(out[[1]], cs(x, 2))

    out <- rowBlockApply(x, rs, mult=2)
    expect_identical(out[[1]], rs(x, 2))
})

test_that("apply on non-pristine DelayedMatrices respects grid construction", {
    # Need to make this non-pristine to avoid fallback to the seed.
    x <- DelayedArray(matrix(runif(10000), ncol=10)) * 2

    # Works with parallelization.
    BPPARAM <- SnowParam(2)

    out <- colBlockApply(x, cs, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_identical(unlist(out), cs(x))

    out <- rowBlockApply(x, rs, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_identical(unlist(out), rs(x))

    # Parallelization and grid size play nice.
    setAutoBlockSize(nrow(x) * 8)

    out <- colBlockApply(x, cs, BPPARAM=BPPARAM)
    expect_identical(length(out), ncol(x))
    expect_identical(unlist(out), cs(x))

    setAutoBlockSize(ncol(x) * 8 * 10)

    out <- rowBlockApply(x, rs, BPPARAM=BPPARAM)
    expect_identical(length(out), as.integer(nrow(x) / 10L))
    expect_identical(unlist(out), rs(x))

    setAutoBlockSize()
})

test_that("apply preserves sparsity in sparse DelayedMatrices", {
    # Need to make this non-pristine to avoid fallback to the seed.
    x <- DelayedArray(Matrix::rsparsematrix(100, 50, density=0.1)) * 2

    # Only one matrix emitted.
    out <- colBlockApply(x, identity, coerce.sparse=FALSE)
    expect_identical(length(out), 1L)
    expect_s4_class(out[[1]], "SVT_SparseArray")

    out <- rowBlockApply(x, identity, coerce.sparse=FALSE)
    expect_identical(length(out), 1L)
    expect_s4_class(out[[1]], "SVT_SparseArray")

    # Coercion to a dgCMatrix works as expected.
    out2 <- colBlockApply(x, identity)
    expect_identical(out2, lapply(out, as, Class="dgCMatrix"))

    sas <- DelayedArray(COO_SparseArray(nzcoo=cbind(1:10, 1:10), 1:10, dim=c(10, 10)))
    out2 <- colBlockApply(sas, identity)
    expect_identical(out2, list(as(seed(sas), "dgCMatrix")))

    # Works with multiple threads. 
    BPPARAM <- SnowParam(2)

    out <- colBlockApply(x, identity, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_true(all(vapply(out, is, class="dgCMatrix", FUN.VALUE=TRUE)))

    out <- rowBlockApply(x, identity, BPPARAM=BPPARAM, coerce.sparse=FALSE)
    expect_identical(length(out), 2L)
    expect_true(all(vapply(out, is, class="SVT_SparseArray", FUN.VALUE=TRUE)))
})

test_that("logical aliases work as expected", {
    x <- matrix(runif(10000), ncol=50)
    setAutoBlockSize(nrow(x) * 8 * 10)

    out <- colBlockApply(x, cs, grid=TRUE)
    expect_equal(length(out), ncol(x)/10L)

    out <- colBlockApply(x, cs, grid=FALSE)
    expect_identical(length(out), 1L)

    setAutoBlockSize(ncol(x) * 8 * 10)

    out <- rowBlockApply(x, rs, grid=TRUE)
    expect_equal(length(out), nrow(x)/10L)

    out <- rowBlockApply(x, rs, grid=FALSE)
    expect_identical(length(out), 1L)

    setAutoBlockSize()
})

test_that("currentViewport works correctly", {
    x <- matrix(runif(10000), ncol=10)
    extractor <- function(x) { 
        DelayedArray::currentViewport()
    }

    # Reports the grid size.
    out <- rowBlockApply(x, extractor)
    expect_identical(length(out), 1L)
    expect_equal(dim(out[[1]]), dim(x))
    expect_s4_class(out[[1]], "ArrayViewport")

    setAutoBlockSize(ncol(x) * 8 * 10)
    out <- rowBlockApply(x, extractor, grid=TRUE)
    expect_equal(length(out), nrow(x)/10L)
    expect_equal(start(out[[100]]), c(991L, 1L))
    expect_equal(end(out[[100]]), c(1000L, 10L))
    expect_s4_class(out[[1]], "ArrayViewport")
    setAutoBlockSize()

    dout <- rowBlockApply(DelayedArray(x) * 2, extractor) # remember, non-pristine for proper checks.
    expect_identical(length(dout), 1L)
    expect_s4_class(dout[[1]], "ArrayViewport")

    # Same handling when there are multiple grid elements.
    out <- rowBlockApply(x, extractor, BPPARAM=SnowParam(2))
    expect_s4_class(out[[1]], "ArrayViewport")
    expect_s4_class(out[[2]], "ArrayViewport")
    expect_false(identical(out[[1]], out[[2]]))

    dout <- rowBlockApply(DelayedArray(x), extractor, BPPARAM=SnowParam(2))
    expect_identical(dout, out)
})
