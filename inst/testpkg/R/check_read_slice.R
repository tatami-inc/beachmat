#' @export
#' @importFrom testthat expect_identical
check_read_slice <- function(FUN, ..., mode) {
    check_read_slice_row(FUN(...), mode)
    check_read_slice_col(FUN(...), mode)
}

#' @importFrom testthat expect_identical
check_read_slice_row <- function(test.mat, mode, FUN="get_row_slice") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    rranges <- spawn_row_ordering(nrow(test.mat))
    NCOL <- ncol(test.mat)
    cbounds <- list(full=c(1L, NCOL), left=c(1L, floor(NCOL/2L)), right=c(ceiling(NCOL/2L), NCOL), middle=sort(sample(NCOL, 2)), single=rep(sample(NCOL, 1), 2))

    for (o in rranges) {
        for (b in cbounds) {
            range <- b[1]:b[2]
            expect_identical(ref[o,range,drop=FALSE], .Call(paste0(FUN, "_", mode), test.mat, o, b, PACKAGE="beachtest"))
        }
    }

    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_read_slice_col <- function(test.mat, mode, FUN="get_col_slice") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    cranges <- spawn_col_ordering(ncol(test.mat))
    NROW <- nrow(test.mat)
    rbounds <- list(full=c(1L, NROW), left=c(1L, floor(NROW/2L)), right=c(ceiling(NROW/2L), NROW), middle=sort(sample(NROW, 2)), single=rep(sample(NROW, 1), 2))

    for (o in cranges) {
        for (b in rbounds) {
            range <- b[1]:b[2]
            expect_identical(ref[range,o,drop=FALSE], .Call(paste0(FUN, "_", mode), test.mat, o, b, PACKAGE="beachtest"))
        }
    }

    return(invisible(NULL))
}
