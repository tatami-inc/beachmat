test_that("internal test package CHECKs correctly", {
    expect_error(rcmdcheck::rcmdcheck(system.file("newtests", package="beachmat"), 
        args=c("--no-multiarch"), error_on="error"), NA)
})
