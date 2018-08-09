#' @export
check_write_type <- function(FUN, ..., mode, out.class=NULL) {
    convertible <- c("logical", "numeric", "integer")
    if (mode %in% convertible) {
        alternative <- setdiff(convertible, mode)

        for (alt.mode in alternative) {
            if (alt.mode=="logical") {
                next # C++ numeric/integer->logical goes via int, which is not quite right.
            }

            test.mat <- FUN(...)
            ref <- as.matrix(test.mat)
            storage.mode(ref) <- alt.mode
        
            ordering <- seq_len(nrow(test.mat))
            out <- .Call(paste0("set_row_", mode, "_to_", alt.mode), test.mat, ordering, PACKAGE="beachtest")
            expect_matrix(ref, out, out.class)

            ordering <- seq_len(ncol(test.mat))
            out <- .Call(paste0("set_col_", mode, "_to_", alt.mode), test.mat, ordering, PACKAGE="beachtest")
            expect_matrix(ref, out, out.class)
        }
    }
    return(invisible(NULL))
}
