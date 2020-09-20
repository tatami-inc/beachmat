test_that("internal test package CHECKs correctly", {
    expect_error(devtools::check(system.file("newtests", package="beachmat"), document=FALSE, error_on="error"), NA)
})
