#' @importFrom DelayedArray defaultGrid
setupDelayedMatrix <- function(mat) {
    grid <- defaultGrid(mat)
    return(list(dim(mat), grid@spacings))
}

realizeDelayedMatrixByRow <- function(mat, I) {
    ind <- (I[1] + 1L):I[2]
    return(as.matrix(mat[ind,,drop=FALSE]))
}

realizeDelayedMatrixByCol <- function(mat, J) { 
    ind <- (J[1] + 1L):J[2]
    return(as.matrix(mat[,ind,drop=FALSE]))
}
