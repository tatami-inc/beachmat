#' @export
#' @importFrom testthat expect_identical
check_write_slice <- function(FUN, ..., mode, out.class) {
    check_write_slice_row(FUN(...), mode, out.class)
    check_write_slice_col(FUN(...), mode, out.class)
}

#' @importFrom testthat expect_identical
check_write_slice_row <- function(test.mat, mode, out.class, FUN="set_row_slice") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    rranges <- spawn_row_ordering(nrow(test.mat))
    cbounds <- spawn_col_bounds(ncol(test.mat))

    for (o in rranges) {
        for (b in cbounds) {
            out <- .Call(paste0(FUN, "_", mode), test.mat, o, b, PACKAGE="beachtest")
            expect_identical(as.character(class(out[[1]])), out.class)

            REF <- ref
            REF[] <- get(mode)(1)
            range <- b[1]:b[2]
            REF[seq_along(o),range] <- ref[o,range]
            expect_equal(REF, as.matrix(out[[1]]))

            expect_equal(REF[o,range,drop=FALSE], out[[2]])
        }
    }
    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_write_slice_col <- function(test.mat, mode, out.class, FUN="set_col_slice") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    cranges <- spawn_col_ordering(ncol(test.mat))
    rbounds <- spawn_row_bounds(nrow(test.mat))

    for (o in cranges) {
        for (b in rbounds) {
            out <- .Call(paste0(FUN, "_", mode), test.mat, o, b, PACKAGE="beachtest")
            expect_identical(as.character(class(out[[1]])), out.class)

            REF <- ref
            REF[] <- get(mode)(1)
            range <- b[1]:b[2]
            REF[range,seq_along(o)] <- ref[range,o]
            expect_equal(REF, as.matrix(out[[1]]))
            
            expect_equal(REF[range,o,drop=FALSE], out[[2]])
        }
    }
    return(invisible(NULL))
}
