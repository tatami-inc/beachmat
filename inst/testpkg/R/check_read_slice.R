#' @export
#' @importFrom testthat expect_identical
check_read_slice <- function(FUN, ..., mode) {
    # Checking row access.
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL
    NROW <- nrow(test.mat)
    NCOL <- ncol(test.mat)

    rranges <- list(forward=seq_len(nrow(test.mat)), random=sample(nrow(test.mat)), subset=sample(nrow(test.mat), nrow(test.mat)/2L))
    rranges$reverse <- rev(rranges$forward)
    cbounds <- list(full=c(1L, NCOL), left=c(1L, floor(NCOL/2L)), right=c(ceiling(NCOL/2L), NCOL), middle=sort(sample(NCOL, 2)), single=rep(sample(NCOL, 1), 2))

    for (o in rranges) {
        for (b in cbounds) {
            range <- b[1]:b[2]
            expect_identical(ref[o,range,drop=FALSE], .Call(paste0("get_row_slice_", mode), test.mat, o, b, PACKAGE="beachtest"))
        }
    }

    # Checking column access.
    test.mat <- FUN(...)
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    cranges <- list(forward=seq_len(ncol(test.mat)), random=sample(ncol(test.mat)), subset=sample(ncol(test.mat), ncol(test.mat)/2L))
    cranges$reverse <- rev(cranges$forward)
    rbounds <- list(full=c(1L, NROW), left=c(1L, floor(NROW/2L)), right=c(ceiling(NROW/2L), NROW), middle=sort(sample(NROW, 2)), single=rep(sample(NROW, 1), 2))

    for (o in cranges) {
        for (b in rbounds) {
            range <- b[1]:b[2]
            expect_identical(ref[range,o,drop=FALSE], .Call(paste0("get_col_slice_", mode), test.mat, o, b, PACKAGE="beachtest"))
        }
    }

    return(invisible(NULL))
}
