SPAWN <- function(nr, nc, mode) {
    mat <- Matrix::rsparsematrix(nr, nc, density=0.3)
    if (mode!=1L) {
        if (mode==0L) mat <- mat != 0
        list(
            as.matrix(mat),
            mat,
            as(mat, "SparseArraySeed")
        )
    } else {
        mat <- round(mat)
        output <- list(
            as.matrix(mat),
            as(mat, "SparseArraySeed")
        )
        storage.mode(output[[1]]) <- "integer"
        storage.mode(output[[2]]@nzdata) <- "integer"
        output
    }
}

CONVERT <- function(x, mode) {
    if (mode==0) {
        storage.mode(x) <- "integer" # as logical conversion goes via integer truncation.
        x <- x != 0L
    } else {
        storage.mode(x) <- c("integer", "double")[mode]
    }
    dimnames(x) <- NULL
    x
}

CHECK_IDENTITY <- function(ref, mat, mode) {
    ref <- CONVERT(ref, mode)
    dimnames(ref) <- NULL

    if (mode==0L) {
        mat <- !!mat # due to the fact that logicals are integers, so non-1 values behave weirdly.
    }
    expect_identical(ref, mat)
}
