#' @importFrom DelayedArray netSubsetAndAperm contentIsPristine nseed seed DelayedArray 
#' @importFrom methods is
setupDelayedMatrix <- function(mat) {
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
            attr(net.sub, "dimmap") <- NULL

            # Creating a matrix for beachmat's API to parse. The logic is that 
            # a seed-only class will not have the "DelayedMatrix" class name when you 
            # wrap it, and it is the wrapped version (e.g., HDF5Matrix, RleMatrix)
            # that constitutes the full matrix. Otherwise, if the wrapped version is
            # just a DelayedMatrix, the seed was a full matrix in the first place.
            wrapped <- DelayedArray(cur.seed)
            if (class(wrapped)[1]=="DelayedMatrix") {
                mat <- cur.seed
            } else {
                mat <- wrapped
            }
            return(list(sub=net.sub, trans=is.trans, mat=mat))
        }
    }

    return(list(sub=NULL, trans=NULL, mat=mat))
}
