# library(testthat); library(beachmat); source("test-tatami-stats.R")

set.seed(1234)
test_that("dimwise sums work as expected", {
    mat <- matrix(runif(1000), 50, 20)
    ptr <- initializeCpp(mat)

    rexpected <- rowSums(mat)
    expect_equal(tatami.sums(ptr, row=TRUE, num.threads=1), rexpected)
    expect_equal(tatami.sums(ptr, row=TRUE, num.threads=2), rexpected)

    cexpected <- colSums(mat)
    expect_equal(tatami.sums(ptr, row=FALSE, num.threads=1), cexpected)
    expect_equal(tatami.sums(ptr, row=FALSE, num.threads=2), cexpected)
})

test_that("dimwise ranges work as expected", {
    mat <- matrix(runif(1000), 25, 40)
    ptr <- initializeCpp(mat)

    rmin <- apply(mat, 1, min)
    rmax <- apply(mat, 1, max)
    rexpected <- list(minimum = rmin, maximum = rmax)
    expect_equal(tatami.ranges(ptr, row=TRUE, num.threads=1), rexpected)
    expect_equal(tatami.ranges(ptr, row=TRUE, num.threads=2), rexpected)

    cmin <- apply(mat, 2, min)
    cmax <- apply(mat, 2, max)
    cexpected <- list(minimum = cmin, maximum = cmax)
    expect_equal(tatami.ranges(ptr, row=FALSE, num.threads=1), cexpected)
    expect_equal(tatami.ranges(ptr, row=FALSE, num.threads=2), cexpected)
})

test_that("dimwise medians work as expected with a more interesting dense matrix", {
    mat <- matrix(runif(1000), 25, 40)
    ptr <- initializeCpp(mat)

    rref <- apply(mat, 1, median)
    expect_equal(tatami.row.medians(ptr, 1), rref)
    expect_equal(tatami.row.medians(ptr, 2), rref)
    cref <- apply(mat, 2, median)
    expect_equal(tatami.column.medians(ptr, 1), cref)
    expect_equal(tatami.column.medians(ptr, 2), cref)
})

test_that("dimwise variances work as expected", { 
    x1 <- Matrix::rsparsematrix(1000, 100, 0.1)
    ptr1 <- initializeCpp(x1)

    out <- tatami.variances(ptr1, row=TRUE, num.threads=1)
    expect_equal(out$mean, rowMeans(x1))
    expect_equal(out$variance, apply(x1, 1, var))

    out <- tatami.variances(ptr1, row=FALSE, num.threads=2)
    expect_equal(out$mean, colMeans(x1))
    expect_equal(out$variance, apply(x1, 2, var))

    # Expected results with dense matrices.
    mat <- matrix(runif(1000), 25, 40)
    dptr <- initializeCpp(mat)

    out <- tatami.variances(dptr, row=TRUE, num.threads=2)
    expect_equal(out$mean, rowMeans(mat))
    expect_equal(out$variance, apply(mat, 1, var))

    out <- tatami.variances(dptr, row=FALSE, num.threads=1)
    expect_equal(out$mean, colMeans(mat))
    expect_equal(out$variance, apply(mat, 2, var))
})

test_that("dimwise NaN counting works as expected", {
    x1 <- Matrix::rsparsematrix(1000, 100, 0.1)

    copy <- x1
    copy[1,1] <- NaN
    copy[2,5] <- NA
    copy[10,5] <- NA
    copy[11,3] <- NA
    copy[11,9] <- NA

    ptr <- initializeCpp(copy)
    expect_equal(tatami.row.nan.counts(ptr, 1), rowSums(is.na(copy)))
    expect_equal(tatami.column.nan.counts(ptr, 1), colSums(is.na(copy)))
})

test_that("dimwise sum by groups work as expected", {
    mat <- matrix(runif(1000), 25, 40)
    ptr <- initializeCpp(mat)

    cnum <- 5
    cgroup <- sample(1:cnum, ncol(mat), replace=TRUE)
    rnum <- 4
    rgroup <- sample(1:rnum, nrow(mat), replace=TRUE)

    rref <- matrix(0, nrow(mat), cnum)
    for (i in 1:5) {
        rref[,i] <- rowSums(mat[, cgroup == i, drop=FALSE])
    }
    expect_equal(tatami.sums.by.group(ptr, cgroup, cnum, row=TRUE, num.threads=1), rref)
    expect_equal(tatami.sums.by.group(ptr, cgroup, cnum, row=TRUE, num.threads=2), rref)

    cref <- matrix(0, rnum, ncol(mat))
    for (i in 1:4) {
        cref[i,] <- colSums(mat[rgroup == i,, drop=FALSE])
    }
    expect_equal(tatami.sums.by.group(ptr, rgroup, rnum, row=FALSE, num.threads=1), cref)
    expect_equal(tatami.sums.by.group(ptr, rgroup, rnum, row=FALSE, num.threads=2), cref)
})

test_that("dimwise variance by groups work as expected", {
    mat <- matrix(runif(1000), 25, 40)
    ptr <- initializeCpp(mat)

    cnum <- 5
    cgroup <- sample(cnum, ncol(mat), replace=TRUE)
    rnum <- 4
    rgroup <- sample(rnum, nrow(mat), replace=TRUE)

    rmean <- rvar <- matrix(0, nrow(mat), cnum)
    for (i in 1:5) {
        submat <- mat[, cgroup == i, drop=FALSE]
        rmean[,i] <- rowMeans(submat)
        rvar[,i] <- apply(submat, 1, var)
    }
    rexpected <- list(mean = rmean, variance = rvar)
    expect_equal(tatami.variances.by.group(ptr, cgroup, cnum, row=TRUE, num.threads=1), rexpected)
    expect_equal(tatami.variances.by.group(ptr, cgroup, cnum, row=TRUE, num.threads=2), rexpected)

    cmean <- cvar <- matrix(0, rnum, ncol(mat))
    for (i in 1:4) {
        submat <- mat[rgroup == i,, drop=FALSE]
        cmean[i,] <- colMeans(submat)
        cvar[i,] <- apply(submat, 2, var)
    }
    cexpected <- list(mean = cmean, variance = cvar)
    expect_equal(tatami.variances.by.group(ptr, rgroup, rnum, row=FALSE, num.threads=1), cexpected)
    expect_equal(tatami.variances.by.group(ptr, rgroup, rnum, row=FALSE, num.threads=2), cexpected)
})

test_that("dimwise median by groups work as expected", {
    mat <- matrix(runif(1000), 25, 40)
    ptr <- initializeCpp(mat)

    cnum <- 5
    cgroup <- sample(1:cnum, ncol(mat), replace=TRUE)
    rnum <- 4
    rgroup <- sample(1:rnum, nrow(mat), replace=TRUE)

    rref <- matrix(0, nrow(mat), cnum)
    for (i in 1:5) {
        rref[,i] <- apply(mat[, cgroup == i, drop=FALSE], 1, median)
    }
    expect_equal(tatami.medians.by.group(ptr, cgroup, cnum, row=TRUE, num.threads=1), rref)
    expect_equal(tatami.medians.by.group(ptr, cgroup, cnum, row=TRUE, num.threads=2), rref)

    cref <- matrix(0, rnum, ncol(mat))
    for (i in 1:4) {
        cref[i,] <- apply(mat[rgroup == i,, drop=FALSE], 2, median)
    }
    expect_equal(tatami.medians.by.group(ptr, rgroup, rnum, row=FALSE, num.threads=1), cref)
    expect_equal(tatami.medians.by.group(ptr, rgroup, rnum, row=FALSE, num.threads=2), cref)
})
