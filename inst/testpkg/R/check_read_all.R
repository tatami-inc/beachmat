#' @export
#' @importFrom testthat expect_identical
check_read_all <- function(FUN, ..., mode) {
    # Checking row access.
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL
    NROW <- nrow(test.mat)
    NCOL <- ncol(test.mat)

    rranges <- list(forward=seq_len(NROW), random=sample(NROW), subset=sample(NROW, NROW/2L))
    rranges$reverse <- rev(rranges$forward)
    for (o in rranges) { 
        expect_identical(ref[o,,drop=FALSE], .Call(paste0("get_row_all_", mode), test.mat, o, PACKAGE="beachtest"))
    }

    # Checking column access.
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    cranges <- list(forward=seq_len(NCOL), random=sample(NCOL), subset=sample(NCOL, NCOL/2L))
    cranges$reverse <- rev(cranges$forward)
    for (o in cranges) {
        expect_identical(ref[,o,drop=FALSE], .Call(paste0("get_col_all_", mode), test.mat, o, PACKAGE="beachtest"))
    }

    # Checking any access.
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    for (ro in rranges) {
        for (co in cranges) { 
            expect_identical(ref[ro,co,drop=FALSE], .Call(paste0("get_single_all_", mode), test.mat, ro, co, PACKAGE="beachtest"))
        }
    }

    return(invisible(NULL))
}
