#' @export
#' @importFrom testthat expect_identical
check_write_varslice <- function(FUN, ..., mode, out.class) {
    check_write_varslice_row(FUN(...), mode, out.class)
    check_write_varslice_col(FUN(...), mode, out.class)
}

#' @importFrom testthat expect_identical
check_write_varslice_row <- function(test.mat, mode, out.class, FUN="set_row_varslice") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL
    NROW <- nrow(ref)	
    NCOL <- ncol(ref)	

    rranges <- spawn_row_ordering(nrow(test.mat))
    cbounds <- spawn_col_bounds(ncol(test.mat))

    for (o in rranges) {
        nentries <- length(o)
        bound1 <- sample(NCOL, nentries, replace=TRUE)
        bound2 <- sample(NCOL, nentries, replace=TRUE)
        cbounds <- cbind(pmin(bound1, bound2), pmax(bound1, bound2))

        out <- .Call(paste0(FUN, "_", mode), test.mat, o, cbounds, PACKAGE="beachtest")
        expect_identical(as.character(class(out[[1]])), out.class)

        REF <- ref
        REF[] <- get(mode)(1)
        for (i in seq_along(o)) {
            range <- cbounds[i,1]:cbounds[i,2]
            REF[i,range] <- ref[o[i], range]
        }
        expect_equal(REF, as.matrix(out[[1]]))

        ref.list <- get_reference_varslice(REF, o, cbounds)
        expect_identical(ref.list, out[[2]])
    }
    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_write_varslice_col <- function(test.mat, mode, out.class, FUN="set_col_varslice") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL
    NROW <- nrow(ref)	
    NCOL <- ncol(ref)	

    cranges <- spawn_col_ordering(ncol(test.mat))
    rbounds <- spawn_row_bounds(nrow(test.mat))

    for (o in cranges) {
        nentries <- length(o)
        bound1 <- sample(NROW, nentries, replace=TRUE)
        bound2 <- sample(NROW, nentries, replace=TRUE)
        rbounds <- cbind(pmin(bound1, bound2), pmax(bound1, bound2))

        out <- .Call(paste0(FUN, "_", mode), test.mat, o, rbounds, PACKAGE="beachtest")
        expect_identical(as.character(class(out[[1]])), out.class)

        REF <- ref
        REF[] <- get(mode)(1)
        for (i in seq_along(o)) {
            range <- rbounds[i,1]:rbounds[i,2]
            REF[range,i] <- ref[range, o[i]]
        }
        expect_equal(REF, as.matrix(out[[1]]))

        ref.list <- get_reference_varslice(REF, o, rbounds, byrow=FALSE)
        expect_identical(ref.list, out[[2]])
    }
    return(invisible(NULL))
}
