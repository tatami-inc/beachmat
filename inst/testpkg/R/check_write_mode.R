#' @export
#' @importFrom testthat expect_identical
check_write_mode <- function(test.mat, expected, simplify=TRUE, preserve.zeroes=TRUE) {
    expect_identical(expected, .Call("set_class_by_sexp", test.mat, simplify, preserve.zeroes, PACKAGE="beachtest"))
}
