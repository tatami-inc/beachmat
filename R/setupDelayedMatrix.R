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

#' @importFrom DelayedArray netSubsetAndAperm contentIsPristine nseed seed DelayedArray 
#' @importFrom methods is
parseDelayedOps <- function(mat) {
    # Finding content-altering modifications.
    no.mod <- contentIsPristine(mat) && nseed(mat)==1L

    if (no.mod) {
        # Finding the net subsetting operations.
        net.sub <- netSubsetAndAperm(mat)

        cur.seed <- seed(mat)
        if (length(dim(cur.seed))==2L) {
            # Determining whether the matrix is transposed (obviously, "dim(mat)"
            # should also be of length 2, for a DelayedMatrix, so 2:1 and NULL
            # are the only possibilities for dimmap).
            dimmap <- attr(net.sub, "dimmap")
            is.trans <- identical(dimmap, 2:1)

            # Creating a matrix for beachmat's API to parse.
            if (is(cur.seed, "HDF5ArraySeed") || is(cur.seed, "RleArraySeed")) {
                mat <- DelayedArray(cur.seed)
            } else {
                mat <- cur.seed
            }
            return(list(sub=net.sub, trans=is.trans, mat=mat))
        }
    }

    return(list(sub=NULL, trans=NULL, mat=mat))
}
