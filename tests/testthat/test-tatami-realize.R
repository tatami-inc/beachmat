# library(testthat); library(beachmat); source("test-tatami-realize.R")

set.seed(9999)
test_that("simple realization works as expected", {
    x1 <- Matrix::rsparsematrix(1000, 100, 0.1)

    d1 <- as.matrix(x1)
    dptr1 <- initializeCpp(d1)
    expect_identical(tatami.realize(dptr1, num.threads=1), d1)
    expect_identical(tatami.realize(dptr1, num.threads=2), d1)

    # To a dgCMatrix.
    ptr1 <- initializeCpp(x1)
    expect_identical(tatami.realize(ptr1, num.threads=1, as.svt=FALSE), x1)
    expect_identical(tatami.realize(ptr1, num.threads=2, as.svt=FALSE), x1)
})

test_that("SVT realization works with column-major sparse matrices", {
    x1 <- Matrix::rsparsematrix(1000, 100, 0.1)

    ptr1 <- initializeCpp(x1)
    svt1 <- tatami.realize(ptr1, num.threads=1, as.svt=TRUE)
    expect_s4_class(svt1, "SVT_SparseMatrix")
    expect_identical(as.matrix(svt1), as.matrix(x1))
    expect_identical(tatami.realize(ptr1, num.threads=2, as.svt=TRUE), svt1)

    # Zeros in the structural non-zeros.    
    haszero <- trunc(x1)
    hptr1 <- initializeCpp(haszero)
    hsvt1 <- tatami.realize(hptr1, num.threads=1, as.svt=TRUE)
    expect_s4_class(hsvt1, "SVT_SparseMatrix")
    expect_identical(as.matrix(hsvt1), as.matrix(haszero))
    expect_identical(tatami.realize(hptr1, num.threads=2, as.svt=TRUE), hsvt1)

    # All 1's. 
    allones <- abs(sign(x1))
    aptr1 <- initializeCpp(allones)
    asvt1 <- tatami.realize(aptr1, num.threads=1, as.svt=TRUE)
    expect_s4_class(asvt1, "SVT_SparseMatrix")
    expect_identical(as.matrix(asvt1), as.matrix(allones))
    expect_identical(tatami.realize(aptr1, num.threads=2, as.svt=TRUE), asvt1)

    # Mix of 1 and zeros. 
    oneorzero <- abs(sign(trunc(x1)))
    optr1 <- initializeCpp(oneorzero)
    osvt1 <- tatami.realize(optr1, num.threads=1, as.svt=TRUE)
    expect_s4_class(osvt1, "SVT_SparseMatrix")
    expect_identical(as.matrix(osvt1), as.matrix(oneorzero))
    expect_identical(tatami.realize(optr1, num.threads=2, as.svt=TRUE), osvt1)

    # Near-empty SVT.    
    nearsparse <- x1
    nearsparse[,sample(ncol(nearsparse), ncol(nearsparse)/2)] <- 0
    nearptr <- initializeCpp(nearsparse)
    nearsvt <- tatami.realize(nearptr, num.threads=1, as.svt=TRUE)
    expect_s4_class(nearsvt, "SVT_SparseMatrix")
    expect_identical(as.matrix(nearsvt), as.matrix(nearsparse))

    # Completely empty SVT.    
    emptysparse <- Matrix::rsparsematrix(20, 5, 0.)
    emptyptr <- initializeCpp(emptysparse)
    emptysvt <- tatami.realize(emptyptr, num.threads=1, as.svt=TRUE)
    expect_s4_class(emptysvt, "SVT_SparseMatrix")
    expect_identical(as.matrix(emptysvt), as.matrix(emptysparse))

    # More columns than rows, just to make sure our buffers are allocated correctly.
    x2 <- Matrix::rsparsematrix(100, 1000, 0.1)

    ptr2 <- initializeCpp(x2)
    svt2 <- tatami.realize(ptr2, num.threads=1, as.svt=TRUE)
    expect_s4_class(svt2, "SVT_SparseMatrix")
    expect_identical(as.matrix(svt2), as.matrix(x2))
    expect_identical(tatami.realize(ptr2, num.threads=2, as.svt=TRUE), svt2)
})

library(DelayedArray)
test_that("SVT realization works with row-major sparse matrices", {
    x1 <- t(DelayedArray(Matrix::rsparsematrix(1000, 100, 0.1)))

    ptr1 <- initializeCpp(x1)
    svt1 <- tatami.realize(ptr1, num.threads=1, as.svt=TRUE)
    expect_s4_class(svt1, "SVT_SparseMatrix")
    expect_identical(as.matrix(svt1), as.matrix(x1))
    expect_identical(tatami.realize(ptr1, num.threads=2, as.svt=TRUE), svt1)

    # Zeros in the structural non-zeros.    
    haszero <- trunc(x1)
    hptr1 <- initializeCpp(haszero)
    hsvt1 <- tatami.realize(hptr1, num.threads=1, as.svt=TRUE)
    expect_s4_class(hsvt1, "SVT_SparseMatrix")
    expect_identical(as.matrix(hsvt1), as.matrix(haszero))
    expect_identical(tatami.realize(hptr1, num.threads=2, as.svt=TRUE), hsvt1)

    # All 1's. 
    allones <- abs(sign(x1))
    aptr1 <- initializeCpp(allones)
    asvt1 <- tatami.realize(aptr1, num.threads=1, as.svt=TRUE)
    expect_s4_class(asvt1, "SVT_SparseMatrix")
    expect_identical(as.matrix(asvt1), as.matrix(allones))
    expect_identical(tatami.realize(aptr1, num.threads=2, as.svt=TRUE), asvt1)

    # Mix of 1 and zeros. 
    oneorzero <- abs(sign(trunc(x1)))
    optr1 <- initializeCpp(oneorzero)
    osvt1 <- tatami.realize(optr1, num.threads=1, as.svt=TRUE)
    expect_s4_class(osvt1, "SVT_SparseMatrix")
    expect_identical(as.matrix(osvt1), as.matrix(oneorzero))
    expect_identical(tatami.realize(optr1, num.threads=2, as.svt=TRUE), osvt1)

    # Near-empty SVT.    
    nearsparse <- Matrix::rsparsematrix(1000, 100, 0.1)
    nearsparse[,sample(ncol(nearsparse), ncol(nearsparse)/2)] <- 0
    nearsparse <- t(DelayedArray(nearsparse))
    nearptr <- initializeCpp(nearsparse)
    nearsvt <- tatami.realize(nearptr, num.threads=1, as.svt=TRUE)
    expect_s4_class(nearsvt, "SVT_SparseMatrix")
    expect_identical(as.matrix(nearsvt), as.matrix(nearsparse))

    # Completely empty SVT.    
    emptysparse <- t(DelayedArray(Matrix::rsparsematrix(20, 50, 0.)))
    emptyptr <- initializeCpp(emptysparse)
    emptysvt <- tatami.realize(emptyptr, num.threads=1, as.svt=TRUE)
    expect_s4_class(emptysvt, "SVT_SparseMatrix")
    expect_identical(as.matrix(emptysvt), as.matrix(emptysparse))

    # More rows than columns, just to make sure our buffers are allocated correctly.
    x2 <- t(DelayedArray(Matrix::rsparsematrix(100, 1000, 0.1)))

    ptr2 <- initializeCpp(x2)
    svt2 <- tatami.realize(ptr2, num.threads=1, as.svt=TRUE)
    expect_s4_class(svt2, "SVT_SparseMatrix")
    expect_identical(as.matrix(svt2), as.matrix(x2))
    expect_identical(tatami.realize(ptr2, num.threads=2, as.svt=TRUE), svt2)
})
