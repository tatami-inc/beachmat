#' @export
#' @import Matrix
setMethod("initializeCpp", "dgCMatrix", function(x, ...) initialize_sparse_matrix(x@x, x@i, x@p, nrow(x), ncol(x), byrow=FALSE))

#' @export
setMethod("initializeCpp", "dgRMatrix", function(x, ...) initialize_sparse_matrix(x@x, x@j, x@p, nrow(x), ncol(x), byrow=TRUE))

####################################################################################
####################################################################################

#' @export
#' @import DelayedArray 
setMethod("initializeCpp", "DelayedMatrix", function(x, ...) {
    initialize(x@seed, ...)
})

#' @export
#' @import DelayedArray 
setMethod("initializeCpp", "DelayedAbind", function(x, ...) {
    collected <- lapply(x@seeds, initialize, ...)
    apply_bind(collected, x@along == 1L)
})

#' @export
setMethod("initializeCpp", "DelayedAperm", function(x, ...) {
    seed <- initialize(x@seed, ...)
    if (x@perm[1] == 2L && x@perm[2] == 1L) {
        apply_transpose(seed)
    } else {
        seed
    }
})

#' @export
setMethod("initializeCpp", "DelayedSubset", function(x, ...) {
    seed <- initialize(x@seed, ...)

    for (i in seq_along(x@index)) {
        idx <- x@index[[i]]
        if (!is.null(idx)) {
            seed <- apply_subset(seed, idx, i == 1L)
        }
    }

    seed
})

#' @export
setMethod("initializeCpp", "DelayedSetDimnames", function(x, ...) {
    initialize(x@seed, ...)
})

####################################################################################
####################################################################################

supported.Ops <- c("+", "-", "*", "/")

.apply_arithmetic <- function(seed, op, val, right, row) {
    if (op == "+") {
        return(apply_addition(seed, val, row))
    } else if (op == "*") {
        return(apply_multiplication(seed, val, row))
    } else if (op == "/") {
        return(apply_division(seed, val, right, row))
    } else if (op == "-") {
        return(apply_subtraction(seed, val, right, row))
    }
}

#' @export
setMethod("initializeCpp", "DelayedUnaryIsoOpWithArgs", function(x, ...) {
    seed <- initialize(x@seed, ...)

    # Figuring out the identity of the operation.
    chosen <- NULL
    for (p in supported.Ops) {
        if (identical(x@OP, get(p, envir=baseenv()))) {
            chosen <- p
            break
        }
    }
    if (is.null(chosen)) {
        stop("unknown operation in ", class(x))
    }

    # Saving the left and right args. There should only be one or the other.
    # as the presence of both is not commutative.
    if (length(x@Rargs) + length(x@Largs) !=1) {
        stop("'DelayedUnaryIsoApWithArgs' should operate on exactly one argument")
    }

    right <- length(x@Rargs) > 0
    if (right) {
        args <- x@Rargs[[1]]
        along <- x@Ralong[1]
    } else {
        args <- x@Largs[[1]]
        along <- x@Lalong[1]
    }
    row <- along == 1L

    .apply_arithmetic(seed, chosen, args, right, row)
})

####################################################################################
####################################################################################

.unary_Math <- function(seed, OP) {
    envir <- environment(OP)
    generic <- envir$`.Generic`

    if (!is.null(generic)) {
        if (generic == "abs") {
            return(apply_abs(seed))
        } else if (generic == "sqrt") {
            return(apply_sqrt(seed))
        } else if (generic == "exp") {
            return(apply_exp(seed))
        } else if (generic == "log1p") {
            return(apply_log1p(seed))
        } else if (generic == "round") {
            if (envir$digits != 0) {
                stop("only 'digits = 0' support for delayed 'round' operations")
            }
            return(apply_round(seed))
        }

        log.base.support <- c(log2=2, log10=10)
        if (generic %in% names(log.base.support)) {
            return(apply_log(seed, log.base.support[[generic]]))
        }
    }

    # Special case for the general case log.
    if (isTRUE(all.equal(as.character(body(OP)), c("log", "a", "base")))) {
        return(apply_log(seed, envir$base))
    }
}

.unary_Ops <- function(seed, OP) {
    envir <- environment(OP)
    generic <- envir$`.Generic`

    if (!is.null(generic)) {
        if (generic %in% supported.Ops) {
            e1 <- envir$e1
            e2 <- envir$e2

            if (missing(e2)) {
                if (generic == "+") {
                    return(seed)
                } else if (generic == "-") {
                    return(apply_multiplication(seed, -1, TRUE))
                } else {
                    stop("second argument can only be missing for unary '+' or '-'")
                }
            } else {
                right <- is(e1, "DelayedArray") # i.e., is the operation applied to the right of the seed?
                val <- if (right) e2 else e1

                # Just guessing that it's applied along the rows, if it's not a scalar.
                return(.apply_arithmetic(seed, generic, val, right, TRUE))
            }
        }
    }
}

#' @export
setMethod("initializeCpp", "DelayedUnaryIsoOpStack", function(x, ...) {
    seed <- initialize(x@seed, ...)
    for (i in seq_along(x@OPS)) { 
        OP <- x@OPS[[i]]
        status <- FALSE 

        if (!status) {
            info <- .unary_Ops(seed, OP)
            if (status <- !is.null(info)) {
                seed <- info
            }
        } 

        if (!status) {
            info <- .unary_Math(seed, OP)
            if (status <- !is.null(info)) {
                seed <- info
            }
        } 

        if (!status) {
            stop("unknown OPS[[", i, "]] function in ", class(x))
        }
    }

    seed
})

