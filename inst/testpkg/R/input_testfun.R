# Creating functions to check each set of matrix inputs.

###############################

#' @importFrom testthat expect_identical
.check_mat <- function(FUN, ..., cxxfun) {
    for (it in seq_len(3)) {
        test.mat <- FUN(...)

        # Specifying the order of the requests.
        if (it==1L) {
            ranges <- list(forward=seq_len(ncol(test.mat))-1L,
                           random=sample(ncol(test.mat))-1L)
            ranges$reverse <- rev(ranges$forward)
        } else if (it==2L) {
            ranges <- list(forward=seq_len(nrow(test.mat))-1L,
                           random=sample(nrow(test.mat))-1L)
            ranges$reverse <- rev(ranges$forward)
        } else {
            ranges <- list(integer(0)) # doesn't matter.
        }

        # Checking that re-ordered requests behave as expected.
        ref <- as.matrix(test.mat)
        dimnames(ref) <- NULL
        for (ordering in ranges) {
            if (it==1L) {
                ref2 <- ref[,ordering+1L,drop=FALSE]
            } else if(it==2L) {
                ref2 <- ref[ordering+1L,,drop=FALSE]
            } else {
                ref2 <- ref
            }
            expect_identical(ref2, .Call(cxxfun, test.mat, it, ordering))
        }
    }
    return(invisible(NULL))
}

#' @export
check_integer_mat <- function(FUN, ...) {
    .check_mat(FUN=FUN, ..., cxxfun=cxx_test_integer_access)
}

#' @export
check_character_mat <- function(FUN, ...) {
    .check_mat(FUN=FUN, ..., cxxfun=cxx_test_character_access)
}

#' @export
check_numeric_mat <- function(FUN, ...) {
    .check_mat(FUN=FUN, ..., cxxfun=cxx_test_numeric_access)
}

#' @export
check_logical_mat <- function(FUN, ...) {
    .check_mat(FUN=FUN, ..., cxxfun=cxx_test_logical_access)
}

###############################

#' @importFrom testthat expect_identical
.check_slices <- function(FUN, ..., by.row, by.col, cxxfun) {
    for (x in by.row) {
        rx <- range(x)

        for (y in by.col) {
            ry <- range(y)
            test.mat <- FUN(...) 
            ref <- as.matrix(test.mat[x, y, drop=FALSE])
            dimnames(ref) <- NULL

            for (i in 1:2) {
                expect_identical(ref, .Call(cxxfun, test.mat, i, rx, ry))
            }
        }
    }
    return(invisible(NULL))
}

#' @export
check_integer_slice <- function(FUN, ..., by.row, by.col) {
    .check_slices(FUN=FUN, ..., by.row=by.row, by.col=by.col, cxxfun=cxx_test_integer_slice)
}

#' @export
check_character_slice <- function(FUN, ..., by.row, by.col) {
    .check_slices(FUN=FUN, ..., by.row=by.row, by.col=by.col, cxxfun=cxx_test_character_slice)
}

#' @export
check_numeric_slice <- function(FUN, ..., by.row, by.col) {
    .check_slices(FUN=FUN, ..., by.row=by.row, by.col=by.col, cxxfun=cxx_test_numeric_slice)
}

#' @export
check_logical_slice <- function(FUN, ..., by.row, by.col) {
    .check_slices(FUN=FUN, ..., by.row=by.row, by.col=by.col, cxxfun=cxx_test_logical_slice)
}

###############################

#' @importFrom testthat expect_identical
.check_const_mat <- function(FUN, ..., cxxfun) {
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL
    expect_identical(ref, .Call(cxxfun, test.mat))
    return(invisible(NULL))
}

#' @export
check_integer_const_mat <- function(FUN, ...) {
    .check_const_mat(FUN=FUN, ..., cxxfun=cxx_test_integer_const_access)
}

#' @export
check_character_const_mat <- function(FUN, ...) {
    .check_const_mat(FUN=FUN, ..., cxxfun=cxx_test_character_const_access)
}

#' @export
check_numeric_const_mat <- function(FUN, ...) {
    .check_const_mat(FUN=FUN, ..., cxxfun=cxx_test_numeric_const_access)
}

#' @export
check_logical_const_mat <- function(FUN, ...) {
    .check_const_mat(FUN=FUN, ..., cxxfun=cxx_test_logical_const_access)
}

#' @importFrom testthat expect_identical
.check_const_slices <- function(FUN, ..., by.row, cxxfun) {
    for (x in by.row) {
        rx <- range(x)
        test.mat <- FUN(...)
        ref <- as.matrix(test.mat)
        dimnames(ref) <- NULL
        expect_identical(ref[x,,drop=FALSE], .Call(cxxfun, test.mat, rx))
    }
    return(invisible(NULL))
}

#' @export
check_integer_const_slice <- function(FUN, ..., by.row) {
    .check_const_slices(FUN=FUN, ..., by.row=by.row, cxxfun=cxx_test_integer_const_slice)
}

#' @export
check_character_const_slice <- function(FUN, ..., by.row) {
    .check_const_slices(FUN=FUN, ..., by.row=by.row, cxxfun=cxx_test_character_const_slice)
}

#' @export
check_numeric_const_slice <- function(FUN, ..., by.row) {
    .check_const_slices(FUN=FUN, ..., by.row=by.row, cxxfun=cxx_test_numeric_const_slice)
}

#' @export
check_logical_const_slice <- function(FUN, ..., by.row) {
    .check_const_slices(FUN=FUN, ..., by.row=by.row, cxxfun=cxx_test_logical_const_slice)
}

###############################

#' @importFrom testthat expect_identical
.check_indexed_mat <- function(FUN, ..., cxxfun) {
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL
    testthat::expect_identical(ref, .Call(cxxfun, test.mat))
    return(invisible(NULL))
}

#' @export
check_integer_indexed_mat <- function(FUN, ...) {
    .check_indexed_mat(FUN=FUN, ..., cxxfun=cxx_test_integer_indexed_access)
}

#' @export
check_numeric_indexed_mat <- function(FUN, ...) {
    .check_indexed_mat(FUN=FUN, ..., cxxfun=cxx_test_numeric_indexed_access)
}

#' @export
check_logical_indexed_mat <- function(FUN, ...) {
    .check_indexed_mat(FUN=FUN, ..., cxxfun=cxx_test_logical_indexed_access)
}

#' @export
check_character_indexed_mat <- function(FUN, ...) {
    .check_indexed_mat(FUN=FUN, ..., cxxfun=cxx_test_character_indexed_access)
}

#' @importFrom testthat expect_identical
.check_indexed_slices <- function(FUN, ..., by.row, cxxfun) {
    for (x in by.row) {
        rx <- range(x)

        test.mat <- FUN(...) 
        ref <- as.matrix(test.mat[x,, drop=FALSE])
        dimnames(ref) <- NULL
        expect_identical(ref, .Call(cxxfun, test.mat, rx))
    }
    return(invisible(NULL))
}

#' @export
check_integer_indexed_slice <- function(FUN, ..., by.row) {
    .check_indexed_slices(FUN=FUN, ..., by.row=by.row, cxxfun=cxx_test_integer_indexed_slice)
}

#' @export
check_numeric_indexed_slice <- function(FUN, ..., by.row) {
    .check_indexed_slices(FUN=FUN, ..., by.row=by.row, cxxfun=cxx_test_numeric_indexed_slice)
}

#' @export
check_logical_indexed_slice <- function(FUN, ..., by.row) {
    .check_indexed_slices(FUN=FUN, ..., by.row=by.row, cxxfun=cxx_test_logical_indexed_slice)
}

#' @export
check_character_indexed_slice <- function(FUN, ..., by.row) {
    .check_indexed_slices(FUN=FUN, ..., by.row=by.row, cxxfun=cxx_test_character_indexed_slice)
}

###############################

#' @export
#' @importFrom testthat expect_identical
check_type <- function(FUN, ..., expected) {
    testthat::expect_identical(expected, .Call(cxx_test_type_check, FUN(...)))
}

#' @importFrom testthat expect_identical
.check_conversion <- function(FUN, ..., cxxfun, rfun) {
    for (it in 1:2) {
        test.mat <- FUN(...)
        ref <- as.matrix(test.mat)
        dimnames(ref) <- NULL
        testthat::expect_identical(rfun(ref), .Call(cxxfun, test.mat, it))
    }
    return(invisible(NULL))
}

#' @export
check_numeric_conversion  <- function(FUN, ...) {
    .check_conversion(FUN=FUN, ..., cxxfun=cxx_test_numeric_to_integer,
                      rfun=function(x) {
                          storage.mode(x) <- "integer" 
                          return(x)
                      })
}

#' @export
check_integer_conversion <- function(FUN, ...) {
    .check_conversion(FUN=FUN, ..., cxxfun=cxx_test_integer_to_numeric, 
                      rfun=function(x) { 
                          storage.mode(x) <- "double"
                          return(x)
                      })
}

#' @export
check_logical_conversion <- function(FUN, ...) {
    .check_conversion(FUN=FUN, ..., cxxfun=cxx_test_logical_to_numeric, 
                      rfun=function(x) { 
                          storage.mode(x) <- "double"
                          return(x)
                      })
    .check_conversion(FUN=FUN, ..., cxxfun=cxx_test_logical_to_integer, 
                      rfun=function(x) { 
                          storage.mode(x) <- "integer"
                          return(x)
                      })
}

###############################

#' @importFrom testthat expect_true expect_error
.check_edge_errors <- function(x, cxxfun) {
    expect_true(.Call(cxxfun, x, 0L))
    expect_error(.Call(cxxfun, x, 1L), "row index out of range")
    expect_error(.Call(cxxfun, x, -1L), "column index out of range")
    expect_error(.Call(cxxfun, x, 2L), "column start index is greater than column end index")
    expect_error(.Call(cxxfun, x, -2L), "row start index is greater than row end index")
    expect_error(.Call(cxxfun, x, 3L), "column end index out of range")
    expect_error(.Call(cxxfun, x, -3L), "row end index out of range")
    return(invisible(NULL))
}

#' @export
check_integer_edge_errors <- function(FUN, ...) {
    .check_edge_errors(FUN(...), cxxfun=cxx_test_integer_edge)
}

#' @export
check_logical_edge_errors <- function(FUN, ...) {
    .check_edge_errors(FUN(...), cxxfun=cxx_test_logical_edge)
}

#' @export
check_numeric_edge_errors <- function(FUN, ...) {
    .check_edge_errors(FUN(...), cxxfun=cxx_test_numeric_edge)
}

#' @export
check_character_edge_errors <- function(FUN, ...) {
    .check_edge_errors(FUN(...), cxxfun=cxx_test_character_edge)
}
