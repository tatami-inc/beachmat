#' @export
#' @importFrom testthat expect_true expect_error expect_identical
check_write_errors <- function(FUN, ..., mode, out.class) {
    x <- FUN(...)
    cxxfun <- paste0("set_errors_", mode)

    for (reget in c(FALSE, TRUE)) {
        out <- .Call(cxxfun, x, 0L, reget, PACKAGE="beachtest")
        expect_identical(as.character(class(out)), out.class)
        expect_true(all(as.matrix(out==get(mode)(1))))
    
        expect_error(.Call(cxxfun, x, 1L, reget, PACKAGE="beachtest"), "row index out of range")
        expect_error(.Call(cxxfun, x, -1L, reget, PACKAGE="beachtest"), "column index out of range")
        expect_error(.Call(cxxfun, x, 2L, reget, PACKAGE="beachtest"), "column start index is greater than column end index")
        expect_error(.Call(cxxfun, x, -2L, reget, PACKAGE="beachtest"), "row start index is greater than row end index")
        expect_error(.Call(cxxfun, x, 3L, reget, PACKAGE="beachtest"), "column end index out of range")
        expect_error(.Call(cxxfun, x, -3L, reget, PACKAGE="beachtest"), "row end index out of range")
    }
    return(invisible(NULL))
}


