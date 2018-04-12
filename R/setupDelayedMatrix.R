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

#' @importFrom DelayedArray netSubsetAndAperm contentIsPristine nseed seed 
#' @importFrom methods is
parseDelayedOps <- function(mat) {
    # Finding content-altering modifications.
    no.mod <- contentIsPristine(mat) && nseed(mat)==1L

    net.sub <- is.trans <- NULL
    if (no.mod) {
        # Finding the net subsetting operations.
        net.sub <- netSubsetAndAperm(mat)
        dimmap <- attr(net.sub, "dimmap")

        # Determining whether the matrix is tranposed. 
        if (identical(dimmap, NULL) || identical(dimmap, 1:2)) {
            is.trans <- FALSE
        } else  if (identical(dimmap, 2:1)) {
            is.trans <- TRUE
        } else {
            no.mod <- FALSE # cannot extract directly from the seed matrix.
        }
    }
      
    # Creating a matrix for beachmat's API to parse.
    if (no.mod) {
        mat <- seed(mat)
        if (is(mat, "HDF5ArraySeed") || is(mat, "RleArraySeed")) {
            mat <- DelayedArray(mat)
        }
    }

    return(list(sub=net.sub, trans=is.trans, mat=mat))
}
