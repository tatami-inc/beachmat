#' @export
#' @importFrom testthat expect_identical
check_read_varslice <- function(FUN, ..., mode) {
    # Checking row access.
    test.mat <- FUN(...)
    NROW <- nrow(test.mat)
    NCOL <- ncol(test.mat)

    rranges <- list(forward=seq_len(NROW), random=sample(NROW), subset=sample(NROW, NROW/2L))
    rranges$reverse <- rev(rranges$forward)

    for (o in rranges) {
        nentries <- length(o)
        bound1 <- sample(NCOL, nentries, replace=TRUE)
        bound2 <- sample(NCOL, nentries, replace=TRUE)
        cbounds <- cbind(pmin(bound1, bound2), pmax(bound1, bound2))

        ref <- vector("list", nentries)
        for (i in seq_len(nentries)) {
            ref[[i]] <- test.mat[o[i], cbounds[i,1]:cbounds[i,2]]
        }
        expect_identical(ref, .Call(paste0("get_row_varslice_", mode), test.mat, o, cbounds, PACKAGE="beachtest"))
    }

    # Checking column access.
    test.mat <- FUN(...)
    cranges <- list(forward=seq_len(NCOL), random=sample(NCOL), subset=sample(NCOL, NCOL/2L))
    cranges$reverse <- rev(cranges$forward)

    for (o in cranges) {
        nentries <- length(o)
        bound1 <- sample(NROW, nentries, replace=TRUE)
        bound2 <- sample(NROW, nentries, replace=TRUE)
        rbounds <- cbind(pmin(bound1, bound2), pmax(bound1, bound2))

        ref <- vector("list", nentries)
        for (i in seq_len(nentries)) {
            ref[[i]] <- test.mat[rbounds[i,1]:rbounds[i,2], o[i]]
        }
        expect_identical(ref, .Call(paste0("get_col_varslice_", mode), test.mat, o, rbounds, PACKAGE="beachtest"))
    }

    return(invisible(NULL))
}
