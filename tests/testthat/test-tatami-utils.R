# library(testthat); library(beachmat); source("test-tatami-utils.R")

# Many of these checks are largely redundant with the initializeCpp checks,
# but we'll call them anyway because they are user-exposed.

library(DelayedArray)
set.seed(1000)
x1 <- Matrix::rsparsematrix(1000, 100, 0.1)
x2 <- Matrix::rsparsematrix(1000, 100, 0.1)

test_that("basic methods work as expected", {
    ptr1 <- initializeCpp(x1)
    expect_equal(tatami.dim(ptr1), dim(x1))

    expect_equal(tatami.row(ptr1, 1), x1[1,])
    expect_equal(tatami.row(ptr1, 1000), x1[1000,])
    expect_equal(tatami.column(ptr1, 1), x1[,1])
    expect_equal(tatami.column(ptr1, 100), x1[,100])

    rref <- Matrix::rowSums(x1)
    expect_equal(tatami.row.sums(ptr1, 1), rref)
    expect_equal(tatami.row.sums(ptr1, 2), rref)
    cref <- Matrix::colSums(x1)
    expect_equal(tatami.column.sums(ptr1, 1), cref)
    expect_equal(tatami.column.sums(ptr1, 2), cref)
})

test_that("bind works as expected", {
    ptr1 <- initializeCpp(x1)
    ptr2 <- initializeCpp(x2)

    combined <- tatami.bind(list(ptr1, ptr2), by.row=TRUE)
    expect_equal(tatami.dim(combined), c(2000L, 100L))
    expect_equal(tatami.column(combined, 1), c(x1[,1], x2[,1]))

    combined <- tatami.bind(list(ptr1, ptr2), by.row=FALSE)
    expect_equal(tatami.dim(combined), c(1000L, 200L))
    expect_equal(tatami.row(combined, 1), c(x1[1,], x2[1,]))
})

test_that("transpose works as expected", {
    ptr1 <- initializeCpp(x1)

    transposed <- tatami.transpose(ptr1)
    expect_equal(tatami.dim(transposed), c(100L, 1000L))
    expect_equal(tatami.row(transposed, 1), x1[,1])
})

test_that("subset works as expected", {
    ptr1 <- initializeCpp(x1)
    subbed <- tatami.subset(ptr1, 100:200, by.row=TRUE)
    expect_equal(tatami.dim(subbed), c(101L, 100L))
    expect_equal(tatami.row(subbed, 1), x1[100,])

    subbed <- tatami.subset(ptr1, 50:10, by.row=FALSE)
    expect_equal(tatami.dim(subbed), c(1000L, 41L))
    expect_equal(tatami.row(subbed, 1), x1[1,50:10])
})

test_that("arithmetic works as expected", {
    ptr1 <- initializeCpp(x1)

    added <- tatami.arith(ptr1, "+", 10.5, by.row=TRUE, right=FALSE)
    expect_equal(tatami.dim(added), dim(x1))
    expect_equal(tatami.row(added, 1), x1[1,] + 10.5)

    added <- tatami.arith(ptr1, "+", 1000:1, by.row=TRUE, right=FALSE)
    expect_equal(tatami.dim(added), dim(x1))
    expect_equal(tatami.row(added, 1), x1[1,] + 1000)

    added <- tatami.arith(ptr1, "+", 100:1, by.row=FALSE, right=FALSE)
    expect_equal(tatami.dim(added), dim(x1))
    expect_equal(tatami.row(added, 1), x1[1,] + 100:1)

    div <- tatami.arith(ptr1, "/", 10.5, by.row=TRUE, right=TRUE)
    expect_equal(tatami.dim(div), dim(x1))
    expect_equal(tatami.row(div, 1), x1[1,] / 10.5)

    div <- tatami.arith(ptr1, "/", 1000:1, by.row=TRUE, right=TRUE)
    expect_equal(tatami.dim(div), dim(x1))
    expect_equal(tatami.row(div, 1), x1[1,] / 1000)

    div <- tatami.arith(ptr1, "/", 100:1, by.row=FALSE, right=TRUE)
    expect_equal(tatami.dim(div), dim(x1))
    expect_equal(tatami.row(div, 1), x1[1,] / 100:1)

    div <- tatami.arith(ptr1, "/", 100:1, by.row=FALSE, right=FALSE)
    expect_equal(tatami.dim(div), dim(x1))
    expect_equal(tatami.row(div, 1), 100:1 / x1[1,])

    expect_error(tatami.arith(ptr1, "foo", 100:1, by.row=FALSE, right=FALSE), "unknown operation")
})

test_that("comparison works as expected", {
    ptr1 <- initializeCpp(x1)

    comp <- tatami.compare(ptr1, ">", 10.5, by.row=TRUE, right=TRUE)
    expect_equal(tatami.dim(comp), dim(x1))
    expect_equal(tatami.row(comp, 1), as.double(x1[1,] > 10.5))

    rthreshold <- 1000:1/1000
    comp <- tatami.compare(ptr1, ">", rthreshold, by.row=TRUE, right=TRUE)
    expect_equal(tatami.dim(comp), dim(x1))
    expect_equal(tatami.row(comp, 1), as.double(x1[1,] > rthreshold[1]))

    comp <- tatami.compare(ptr1, ">", rthreshold, by.row=TRUE, right=FALSE)
    expect_equal(tatami.dim(comp), dim(x1))
    expect_equal(tatami.row(comp, 1), as.double(x1[1,] < rthreshold[1]))

    cthreshold <- 1:100/100
    comp <- tatami.compare(ptr1, ">", cthreshold, by.row=FALSE, right=TRUE)
    expect_equal(tatami.dim(comp), dim(x1))
    expect_equal(tatami.row(comp, 1), as.double(x1[1,] > cthreshold))
})

test_that("logic works as expected", {
    ptr1 <- initializeCpp(x1)

    logic <- tatami.logic(ptr1, "&", FALSE, by.row=TRUE)
    expect_equal(tatami.dim(logic), dim(x1))
    expect_equal(tatami.row(logic, 1), as.double(x1[1,] & FALSE))

    rthreshold <- rep(c(TRUE, FALSE), length.out=1000)
    logic <- tatami.logic(ptr1, "&", rthreshold, by.row=TRUE)
    expect_equal(tatami.dim(logic), dim(x1))
    expect_equal(tatami.row(logic, 1), as.double(x1[1,] & rthreshold[1]))

    logic <- tatami.logic(ptr1, "&", rthreshold, by.row=TRUE)
    expect_equal(tatami.dim(logic), dim(x1))
    expect_equal(tatami.row(logic, 1), as.double(x1[1,] & rthreshold[1]))

    cthreshold <- rep(c(TRUE, FALSE), length.out=100)
    logic <- tatami.logic(ptr1, "&", cthreshold, by.row=FALSE)
    expect_equal(tatami.dim(logic), dim(x1))
    expect_equal(tatami.row(logic, 1), as.double(x1[1,] & cthreshold))
})

test_that("rounding works as expected", {
    ptr1 <- initializeCpp(x1)
    rounded <- tatami.round(ptr1)
    expect_equal(tatami.dim(rounded), dim(x1))
    expect_equal(tatami.row(rounded, 1), round(x1[1,]))
})

test_that("logging works as expected", {
    y <- abs(x1) + 1 
    ptr1 <- initializeCpp(y)
    logged <- tatami.log(ptr1, 10)
    expect_equal(tatami.dim(logged), dim(x1))
    expect_equal(tatami.row(logged, 1), log10(y[1,]))
})

test_that("math works as expected", {
    ptr1 <- initializeCpp(x1)
    expped <- tatami.math(ptr1, "exp")
    expect_equal(tatami.dim(expped), dim(x1))
    expect_equal(tatami.row(expped, 1), exp(x1[1,]))
})

test_that("NOT works as expected", {
    ptr1 <- initializeCpp(x1)
    not <- tatami.not(ptr1)
    expect_equal(tatami.dim(not), dim(x1))
    expect_equal(tatami.row(not, 1), as.double(!(x1[1,])))
})

test_that("binary works as expected", {
    ptr1 <- initializeCpp(x1)
    ptr2 <- initializeCpp(x2)
    bin <- tatami.binary(ptr1, ptr2, "+")
    expect_equal(tatami.dim(bin), dim(x1))
    expect_equal(tatami.row(bin, 1), x1[1,] + x2[1,])
})
