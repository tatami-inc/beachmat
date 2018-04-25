# Checks for proper evaluation of delayed operations.
# library(testthat); library(beachmat); source("test-delayed.R")

set.seed(1000)
library(DelayedArray)
delayed_ord <- DelayedArray(matrix(runif(10000), 20, 500))

library(Matrix)
delayed_sparse <- DelayedArray(rsparsematrix(100, 20, density=0.1))

library(HDF5Array)
full_hdf5 <- as(matrix(rnorm(900), 90, 10), "HDF5Array")
delayed_hdf5 <- full_hdf5[1:20,]

full_rle <- as(matrix(rpois(2000, lambda=0.5), 40, 50), "RleArray")
delayed_rle <- full_rle[, 20:11]

test_that("realization methods work correctly", {
    out <- beachmat:::setupDelayedMatrix(delayed_ord)
    expect_identical(out[[1]], dim(delayed_ord))
    expect_identical(length(out[[2]]), 2L)

    out <- beachmat:::setupDelayedMatrix(delayed_sparse)
    expect_identical(out[[1]], dim(delayed_sparse))
    expect_identical(length(out[[2]]), 2L)

    out <- beachmat:::setupDelayedMatrix(delayed_hdf5)
    expect_identical(out[[1]], dim(delayed_hdf5))
    expect_identical(length(out[[2]]), 2L)

    # Realizing by row.
    real <- beachmat:::realizeDelayedMatrixByRow(delayed_ord, c(0, 5))
    expect_equal(real, as.matrix(delayed_ord[1:5,]))

    real <- beachmat:::realizeDelayedMatrixByRow(delayed_sparse, c(1, 6))
    expect_equal(real, as.matrix(delayed_sparse[2:6,]))

    real <- beachmat:::realizeDelayedMatrixByRow(delayed_ord, c(3, 10))
    expect_equal(real, as.matrix(delayed_ord[4:10,]))

    # Realizing by column.
    real <- beachmat:::realizeDelayedMatrixByCol(delayed_ord, c(0, 5))
    expect_equal(real, as.matrix(delayed_ord[,1:5]))

    real <- beachmat:::realizeDelayedMatrixByCol(delayed_sparse, c(1, 6))
    expect_equal(real, as.matrix(delayed_sparse[,2:6]))

    real <- beachmat:::realizeDelayedMatrixByCol(delayed_ord, c(3, 10))
    expect_equal(real, as.matrix(delayed_ord[,4:10]))
})

test_that("delayed operations parsing works correctly", {
    # Checking out all transposing and subsetting.
    parsed <- beachmat:::parseDelayedOps(delayed_ord)
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_identical(class(parsed$mat), "matrix")

    mod1 <- t(delayed_ord)
    parsed <- beachmat:::parseDelayedOps(mod1)
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_identical(class(parsed$mat), "matrix")

    mod2 <- delayed_ord[20:1,]
    parsed <- beachmat:::parseDelayedOps(mod2)
    expect_identical(parsed$sub, list(20:1, NULL))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_identical(class(parsed$mat), "matrix")

    mod3 <- t(mod2)
    parsed <- beachmat:::parseDelayedOps(mod3)
    expect_identical(parsed$sub, list(20:1, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_identical(class(parsed$mat), "matrix")

    mod4 <- mod1[,20:1]
    parsed <- beachmat:::parseDelayedOps(mod4)
    expect_identical(parsed$sub, list(20:1, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_identical(class(parsed$mat), "matrix")

    # Checking out an actual opertion.
    xmod <- delayed_ord + 1
    parsed <- beachmat:::parseDelayedOps(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)

    xmod <- cbind(delayed_ord, delayed_ord)
    parsed <- beachmat:::parseDelayedOps(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)

    xmod <- delayed_ord[1:10,] + 1:10
    parsed <- beachmat:::parseDelayedOps(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)

    xmod <- t(delayed_ord) * 2
    parsed <- beachmat:::parseDelayedOps(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)
 
    # Checking out array inputs with >2 dimensions.
    whee <- DelayedArray(array(runif(6000), c(10, 20, 30)))[,,1]
    parsed <- beachmat:::parseDelayedOps(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)

    # Checking out what happens with the seeds.
    parsed <- beachmat:::parseDelayedOps(delayed_hdf5)
    expect_identical(parsed$sub, list(1:20, NULL))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, full_hdf5)
    expect_s4_class(parsed$mat, "HDF5Matrix")

    parsed <- beachmat:::parseDelayedOps(t(full_hdf5))
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, full_hdf5)
    expect_s4_class(parsed$mat, "HDF5Matrix")

    parsed <- beachmat:::parseDelayedOps(delayed_rle)
    expect_identical(parsed$sub, list(NULL, 20:11))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, full_rle)
    expect_s4_class(parsed$mat, "RleMatrix")

    parsed <- beachmat:::parseDelayedOps(t(full_rle))
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, full_rle)
    expect_s4_class(parsed$mat, "RleMatrix")

    # Checking that sparse entities are not evaluated.
    parsed <- beachmat:::parseDelayedOps(delayed_sparse)
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, seed(delayed_sparse))
    expect_s4_class(parsed$mat, "dgCMatrix")

    parsed <- beachmat:::parseDelayedOps(t(delayed_sparse[9:5,3:8]))
    expect_identical(parsed$sub, list(9:5, 3:8))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, seed(delayed_sparse))
    expect_s4_class(parsed$mat, "dgCMatrix")
})
