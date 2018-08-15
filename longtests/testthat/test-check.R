test_that("internal test package CHECKs correctly", {
    out <- devtools::check(system.file("testpkg", package="beachmat"), document=FALSE)
    expect_identical(out$error, character(0))
})
