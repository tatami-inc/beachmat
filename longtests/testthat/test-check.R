test_that("internal test package CHECKs correctly", {
    expect_error(rcmdcheck::rcmdcheck(system.file("testpkg", package="beachmat"), error_on="error"), NA)
})

test_that("extension test package CHECKs correctly", {
    expect_error(rcmdcheck::rcmdcheck(system.file("extensions", package="beachmat"), document=FALSE, error_on="error"), NA)
})
