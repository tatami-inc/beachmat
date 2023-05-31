# This checks the in-memory caching.
# library(testthat); library(beachmat); source("test-checkMemoryCache.R")

# Mocking up a class with some kind of uniquely identifying aspect.
setClass("UnknownMatrix", slots=c(contents="dgCMatrix", uuid="character"))

# Defining our initialization method.
setMethod("initializeCpp", "UnknownMatrix", function(x, ..., memorize=FALSE) {
    if (memorize) {
        checkMemoryCache("my_package", x@uuid, function() initializeCpp(x@contents))
    } else {
        initializeCpp(x@contents)
    }
})

test_that("in-memory caching and flushing works as expected", {
    X <- new("UnknownMatrix", contents=Matrix::rsparsematrix(10, 10, 0.1), uuid=as.character(sample(1e8, 1)))

    ptr0 <- initializeCpp(X) 
    ptr1 <- initializeCpp(X, memorize=TRUE)
    expect_false(identical(capture.output(print(ptr0)), capture.output(print(ptr1))))

    ptr2 <- initializeCpp(X, memorize=TRUE)
    expect_identical(capture.output(print(ptr1)), capture.output(print(ptr2)))

    flushMemoryCache()
    ptr3 <- initializeCpp(X, memorize=TRUE)
    expect_false(identical(capture.output(print(ptr2)), capture.output(print(ptr3))))
})
