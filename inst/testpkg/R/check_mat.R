#' @importFrom testthat expect_identical
check_mat <- function(FUN, ..., mode) {
    cxx_get_row_all <- get(paste0("cxx_get_row_all_", mode))
    cxx_get_col_all <- get(paste0("cxx_get_col_all_", mode))
    cxx_get_single_all <- get(paste0("cxx_get_single_all_", mode))

    # Checking row access.
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    rranges <- list(forward=seq_len(nrow(test.mat)), random=sample(nrow(test.mat)), subset=sample(nrow(test.mat), nrow(test.mat)/2L))
    rranges$reverse <- rev(rranges$forward)
    for (o in rranges) { 
        expect_identical(ref[o,,drop=FALSE], .Call(cxx_get_row_all, test.mat, o))
    }

    # Checking column access.
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    cranges <- list(forward=seq_len(ncol(test.mat)), random=sample(ncol(test.mat)), subset=sample(ncol(test.mat), ncol(test.mat)/2L))
    cranges$reverse <- rev(cranges$forward)
    for (o in cranges) {
        expect_identical(ref[,o,drop=FALSE], .Call(cxx_get_col_all, test.mat, o))
    }

    # Checking any access.
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    for (ro in rranges) {
        for (co in cranges) { 
            expect_identical(ref[ro,co,drop=FALSE], .Call(cxx_get_single_all, test.mat, ro, co))
        }
    }

    return(invisible(NULL))
}

#' @export
check_integer_mat <- function(FUN, ...) {
    check_mat(FUN=FUN, ..., mode="integer")
}

#' @export
check_character_mat <- function(FUN, ...) {
    check_mat(FUN=FUN, ..., mode="character")
}

#' @export
check_numeric_mat <- function(FUN, ...) {
    check_mat(FUN=FUN, ..., mode="numeric")
}

#' @export
check_logical_mat <- function(FUN, ...) {
    check_mat(FUN=FUN, ..., mode="logical")
}

