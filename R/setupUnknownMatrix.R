#' @importFrom DelayedArray blockGrid
setupUnknownMatrix <- function(mat) {
    grid <- blockGrid(mat)
    list(dim(mat), grid@spacings)
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
