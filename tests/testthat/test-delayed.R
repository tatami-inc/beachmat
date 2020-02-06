# Checks for proper evaluation of delayed operations.
# library(testthat); library(beachmat); source("test-delayed.R")

set.seed(1000)
library(DelayedArray)
delayed_ord <- DelayedArray(matrix(runif(10000), 20, 500))

library(Matrix)
delayed_sparse <- DelayedArray(rsparsematrix(100, 20, density=0.1))

full_rle <- as(matrix(rpois(2000, lambda=0.5), 40, 50), "RleArray")
delayed_rle <- full_rle[, 20:11]

full_rle2 <- as(matrix(rpois(900, lambda=0.5), 90, 10), "RleArray")
delayed_rle2 <- full_rle2[1:20,]

test_that("delayed operations parsing works correctly", {
    # Checking out all transposing and subsetting.
    parsed <- beachmat:::setupDelayedMatrix(delayed_ord)
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_true(is.matrix(parsed$mat))

    mod1 <- t(delayed_ord)
    parsed <- beachmat:::setupDelayedMatrix(mod1)
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_true(is.matrix(parsed$mat))

    mod2 <- delayed_ord[20:1,]
    parsed <- beachmat:::setupDelayedMatrix(mod2)
    expect_identical(parsed$sub, list(20:1, NULL))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_true(is.matrix(parsed$mat))

    mod3 <- t(mod2)
    parsed <- beachmat:::setupDelayedMatrix(mod3)
    expect_identical(parsed$sub, list(20:1, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_true(is.matrix(parsed$mat))

    mod4 <- mod1[,20:1]
    parsed <- beachmat:::setupDelayedMatrix(mod4)
    expect_identical(parsed$sub, list(20:1, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, seed(delayed_ord))
    expect_true(is.matrix(parsed$mat))

    # Checking out an actual opertion.
    xmod <- delayed_ord + 1
    parsed <- beachmat:::setupDelayedMatrix(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)

    xmod <- BiocGenerics::cbind(delayed_ord, delayed_ord)
    parsed <- beachmat:::setupDelayedMatrix(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)

    xmod <- delayed_ord[1:10,] + 1:10
    parsed <- beachmat:::setupDelayedMatrix(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)

    xmod <- t(delayed_ord) * 2
    parsed <- beachmat:::setupDelayedMatrix(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)
 
    # Checking out array inputs with >2 dimensions.
    whee <- DelayedArray(array(runif(6000), c(10, 20, 30)))[,,1]
    parsed <- beachmat:::setupDelayedMatrix(xmod)
    expect_identical(parsed$sub, NULL)
    expect_identical(parsed$trans, NULL)
    expect_identical(parsed$mat, xmod)

    # Checking out what happens with the seeds.
    parsed <- beachmat:::setupDelayedMatrix(delayed_rle)
    expect_identical(parsed$sub, list(NULL, 20:11))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, full_rle)
    expect_s4_class(parsed$mat, "RleMatrix")

    parsed <- beachmat:::setupDelayedMatrix(t(full_rle))
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, full_rle)
    expect_s4_class(parsed$mat, "RleMatrix")

    parsed <- beachmat:::setupDelayedMatrix(delayed_rle2)
    expect_identical(parsed$sub, list(1:20, NULL))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, full_rle2)
    expect_s4_class(parsed$mat, "RleMatrix")

    parsed <- beachmat:::setupDelayedMatrix(t(full_rle2))
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, full_rle2)
    expect_s4_class(parsed$mat, "RleMatrix")

    # Checking that sparse entities are not evaluated.
    parsed <- beachmat:::setupDelayedMatrix(delayed_sparse)
    expect_identical(parsed$sub, list(NULL, NULL))
    expect_identical(parsed$trans, FALSE)
    expect_identical(parsed$mat, seed(delayed_sparse))
    expect_s4_class(parsed$mat, "dgCMatrix")

    parsed <- beachmat:::setupDelayedMatrix(t(delayed_sparse[9:5,3:8]))
    expect_identical(parsed$sub, list(9:5, 3:8))
    expect_identical(parsed$trans, TRUE)
    expect_identical(parsed$mat, seed(delayed_sparse))
    expect_s4_class(parsed$mat, "dgCMatrix")
})
