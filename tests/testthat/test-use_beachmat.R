# Checks for proper execution of use_beachmat.
context("use_beachmat.R")

# Adapted from
# https://github.com/r-lib/usethis/blob/master/tests/testthat/test-use-rcpp.R
test_that("use_beachmat() requires a package", {
    scoped_temporary_project()
    expect_error(beachmat:::use_beachmat(), "not an R package")
})

test_that("use_beachmat() creates files/dirs, edits DESCRIPTION and .gitignore, and adds Makevars[.win]", {
    if (requireNamespace("desc", quietly = TRUE)) {
        pkg <- scoped_temporary_package()
        capture_output(use_beachmat())
        expect_match(desc::desc_get("LinkingTo", pkg), "Rcpp")
        expect_match(desc::desc_get("LinkingTo", pkg), "Rhdf5lib")
        expect_match(desc::desc_get("LinkingTo", pkg), "beachmat")
        expect_match(desc::desc_get("Imports", pkg), "Rcpp")
        expect_match(desc::desc_get("Suggests", pkg), "beachmat")
        expect_match(desc::desc_get("Suggests", pkg), "HDF5Array")
        expect_match(desc::desc_get("Suggests", pkg), "DelayedArray")
        expect_match(desc::desc_get("Suggests", pkg), "Matrix")
        expect_match(desc::desc_get("SystemRequirements", pkg), "C\\+\\+11")
        expect_true(dir.exists(usethis:::proj_path("src")))
        ignores <- readLines(usethis:::proj_path("src", ".gitignore"))
        expect_true(all(c("*.o", "*.so", "*.dll") %in% ignores))
        expect_true(file.exists(usethis:::proj_path("src", "Makevars")))
        expect_true(file.exists(usethis:::proj_path("src", "Makevars.win")))
    }
})
