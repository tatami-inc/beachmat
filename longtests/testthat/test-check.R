test_that("internal test package CHECKs correctly", {
    expect_error(devtools::check(system.file("testpkg", package="beachmat"), document=FALSE, error_on="error"), NA)
})

test_that("extension test package CHECKs correctly", {
    expect_error(devtools::check(system.file("extensions", package="beachmat"), document=FALSE, error_on="error"), NA)
})
