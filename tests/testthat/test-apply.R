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

    out <- colBlockApply(x, cs, BPPARAM=BPPARAM)
    expect_identical(length(out), ncol(x))
    expect_identical(unlist(out), cs(x))

    setAutoBlockSize(ncol(x) * 8 * 10)

    out <- rowBlockApply(x, rs, BPPARAM=BPPARAM)
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

    out <- colBlockApply(x, cs, BPPARAM=BPPARAM)
    expect_identical(length(out), ncol(x))
    expect_identical(unlist(out), cs(x))

    setAutoBlockSize(ncol(x) * 8 * 10)

    out <- rowBlockApply(x, rs, BPPARAM=BPPARAM)
    expect_identical(length(out), as.integer(nrow(x) / 10L))
    expect_identical(unlist(out), rs(x))

    setAutoBlockSize()
})

test_that("apply works with DelayedMatrices", {
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
})

test_that("apply on DelayedMatrices respects grid construction", {
    x <- DelayedArray(matrix(runif(10000), ncol=10))

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
    x <- DelayedArray(Matrix::rsparsematrix(100, 50, density=0.1))

    # Only one matrix emitted.
    out <- colBlockApply(x, identity)
    expect_identical(length(out), 1L)
    expect_true(is(out[[1]], "SparseArraySeed"))

    out <- rowBlockApply(x, identity)
    expect_identical(length(out), 1L)
    expect_true(is(out[[1]], "SparseArraySeed"))

    # Works with multiple matrices.
    BPPARAM <- SnowParam(2)

    out <- colBlockApply(x, identity, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_true(all(vapply(out, is, class="SparseArraySeed", FUN.VALUE=TRUE)))

    out <- rowBlockApply(x, identity, BPPARAM=BPPARAM)
    expect_identical(length(out), 2L)
    expect_true(all(vapply(out, is, class="SparseArraySeed", FUN.VALUE=TRUE)))
})

test_that("native apply reports the grid in attributes", {
    x <- matrix(runif(10000), ncol=10)

    # Reports the grid size.
    out <- rowBlockApply(x, attr, which="from_grid")
    expect_s4_class(out[[1]], "ArrayGrid")

    dout <- rowBlockApply(DelayedArray(x), attr, which="from_grid")
    expect_identical(dout, out)

    # Same handling when there are multiple grid elenehts.
    out <- rowBlockApply(x, attr, which="from_grid", BPPARAM=SnowParam(2))
    expect_s4_class(out[[1]], "ArrayGrid")
    expect_s4_class(out[[2]], "ArrayGrid")

    dout <- rowBlockApply(x, attr, which="from_grid", BPPARAM=SnowParam(2))
    expect_identical(dout, out)
})
