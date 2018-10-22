#' @importFrom DelayedArray rowGrid colGrid
#' @importFrom BiocGenerics dims
setupUnknownMatrix <- function(mat) {
    if (ncol(mat)) {
        byrow <- dims(rowGrid(mat))[,1]
    } else {
        byrow <- nrow(mat)
    }

    if (nrow(mat)) {
        bycol <- dims(colGrid(mat))[,2]
    } else {
        bycol <- ncol(mat)
    }

    list(dim(mat), 
        c(0L, cumsum(byrow)), 
        c(0L, cumsum(bycol))
    )
}

#' @importFrom BiocGenerics t
realizeByRange <- function(mat, i, j, transpose=FALSE) 
# The first element is assumed to encode the start position, 
# while the second element is presumed to encode the length.
# Transposition is useful for effective column access from rows.
{
    I <- i[1] + seq_len(i[2]) # free conversion to 1-based indexing.
    J <- j[1] + seq_len(j[2])
    mat <- mat[I,J,drop=FALSE]

    if (transpose) {
        mat <- t(mat)
    }
    as.matrix(mat)
}

realizeByRangeIndex <- function(mat, i, J) {
    I <- i[1] + seq_len(i[2])
    as.matrix(mat[I,J,drop=FALSE])
}

realizeByIndexRange <- function(mat, I, j) {
    J <- j[1] + seq_len(j[2])
    as.matrix(mat[I,J,drop=FALSE])
}
