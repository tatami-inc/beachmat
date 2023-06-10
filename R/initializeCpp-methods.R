#' @export
setMethod("initializeCpp", "ANY", function(x, ...) initialize_unknown_matrix(x))

#' @export
setMethod("initializeCpp", "matrix", function(x, ...) initialize_dense_matrix(x, nrow(x), ncol(x)))

#' @export
setMethod("initializeCpp", "dgeMatrix", function(x, ...) initialize_dense_matrix(x@x, nrow(x), ncol(x)))

#' @export
setMethod("initializeCpp", "lgeMatrix", function(x, ...) initialize_dense_matrix(x@x, nrow(x), ncol(x)))

####################################################################################
####################################################################################

#' @export
#' @import Matrix
setMethod("initializeCpp", "dgCMatrix", function(x, ...) initialize_sparse_matrix(x@x, x@i, x@p, nrow(x), ncol(x), byrow=FALSE))

#' @export
setMethod("initializeCpp", "dgRMatrix", function(x, ...) initialize_sparse_matrix(x@x, x@j, x@p, nrow(x), ncol(x), byrow=TRUE))

#' @export
#' @import Matrix
setMethod("initializeCpp", "lgCMatrix", function(x, ...) initialize_sparse_matrix(x@x, x@i, x@p, nrow(x), ncol(x), byrow=FALSE))

#' @export
setMethod("initializeCpp", "lgRMatrix", function(x, ...) initialize_sparse_matrix(x@x, x@j, x@p, nrow(x), ncol(x), byrow=TRUE))

####################################################################################
####################################################################################

#' @export
#' @import DelayedArray 
setMethod("initializeCpp", "DelayedMatrix", function(x, ...) {
    initializeCpp(x@seed, ...)
})

#' @export
#' @import DelayedArray 
setMethod("initializeCpp", "DelayedAbind", function(x, ...) {
    collected <- lapply(x@seeds, initializeCpp, ...)
    apply_delayed_bind(collected, x@along == 1L)
})

#' @export
setMethod("initializeCpp", "DelayedAperm", function(x, ...) {
    seed <- initializeCpp(x@seed, ...)
    if (x@perm[1] == 2L && x@perm[2] == 1L) {
        apply_delayed_transpose(seed)
    } else {
        seed
    }
})

#' @export
setMethod("initializeCpp", "DelayedSubset", function(x, ...) {
    seed <- initializeCpp(x@seed, ...)

    for (i in seq_along(x@index)) {
        idx <- x@index[[i]]
        if (!is.null(idx)) {
            seed <- apply_delayed_subset(seed, idx, i == 1L)
        }
    }

    seed
})

#' @export
setMethod("initializeCpp", "DelayedSetDimnames", function(x, ...) {
    initializeCpp(x@seed, ...)
})

####################################################################################
####################################################################################

#' @export
#' @importClassesFrom SparseArray SVT_SparseMatrix
setMethod("initializeCpp", "SVT_SparseMatrix", function(x, ...) {
    initialize_SVT_SparseMatrix(nr=nrow(x), nc=ncol(x), x)
})

####################################################################################
####################################################################################

supported.Arith1 <- c("+", "*")
supported.Arith2 <- c("-", "/")
supported.Compare <- c("==", ">", "<", ">=", "<=", "!=")
supported.Logic <- c("&", "|")
supported.Ops <- c(supported.Arith1, supported.Arith2, supported.Compare, supported.Logic)

reverse.Compare <- c("=="="==", ">"="<", "<"=">", ">="="<=", "<="=">=", "!="="!=")

.apply_delayed_unary_ops <- function(seed, op, val, right, row) {
    if (op %in% supported.Arith1) {
        return(apply_delayed_associative_arithmetic(seed, val, row, op))
    } else if (op %in% supported.Arith2) {
        return(apply_delayed_nonassociative_arithmetic(seed, val, right, row, op))
    } else if (op %in% supported.Compare) {
        if (!right) { # need to flip the operation if the argument is not on the right.
            op <- reverse.Compare[[op]]
        }
        return(apply_delayed_comparison(seed, val, row, op))
    } else if (op %in% supported.Logic) {
        return(apply_delayed_boolean(seed, val, row, op))
    } else {
        return(NULL)
    }
}

#' @export
setMethod("initializeCpp", "DelayedUnaryIsoOpWithArgs", function(x, ...) {
    seed <- initializeCpp(x@seed, ...)

    # Saving the left and right args. There should only be one or the other.
    # as the presence of both is not commutative.
    if (length(x@Rargs) + length(x@Largs) !=1) {
        stop("'", class(x)[1], "' should operate on exactly one argument")
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

    # Figuring out the identity of the operation.
    chosen <- NULL
    for (p in supported.Ops) {
        if (identical(x@OP, get(p, envir=baseenv()))) {
            chosen <- p
            break
        }
    }

    if (is.null(chosen)) {
        warning("unknown operation in '<", class(x)[1], ">@OP', falling back to an unknown matrix")
        return(initialize_unknown_matrix(x))
    } 

    output <- .apply_delayed_unary_ops(seed, chosen, args, right, row)
    stopifnot(!is.null(output))
    output
})

####################################################################################
####################################################################################

.unary_Math <- function(seed, OP) {
    envir <- environment(OP)
    generic <- envir$`.Generic`

    if (is.null(generic)) {
        # Special case for the general log().
        if (isTRUE(all.equal(as.character(body(OP)), c("log", "a", "base")))) {
            return(apply_delayed_log(seed, envir$base))
        }
        return(NULL)
    }

    if (generic == "abs") {
        return(apply_delayed_abs(seed))
    }

    if (generic == "sqrt") {
        return(apply_delayed_sqrt(seed))
    }

    if (generic == "ceiling") {
        return(apply_delayed_ceiling(seed))
    }

    if (generic == "floor") {
        return(apply_delayed_floor(seed))
    }

    if (generic == "trunc") {
        return(apply_delayed_trunc(seed))
    }

    if (generic == "exp") {
        return(apply_delayed_exp(seed))
    }

    if (generic == "expm1") {
        return(apply_delayed_expm1(seed))
    }

    if (generic == "log1p") {
        return(apply_delayed_log1p(seed))
    }

    if (generic == "acos") {
        return(apply_delayed_acos(seed))
    }

    if (generic == "acosh") {
        return(apply_delayed_acosh(seed))
    }

    if (generic == "asin") {
        return(apply_delayed_asin(seed))
    }

    if (generic == "asinh") {
        return(apply_delayed_asinh(seed))
    }

    if (generic == "atan") {
        return(apply_delayed_atan(seed))
    }

    if (generic == "atanh") {
        return(apply_delayed_atanh(seed))
    }

    if (generic == "cos") {
        return(apply_delayed_cos(seed))
    }

    if (generic == "cosh") {
        return(apply_delayed_cosh(seed))
    }

    if (generic == "sin") {
        return(apply_delayed_sin(seed))
    }

    if (generic == "sinh") {
        return(apply_delayed_sinh(seed))
    }

    if (generic == "tan") {
        return(apply_delayed_tan(seed))
    }

    if (generic == "tanh") {
        return(apply_delayed_tanh(seed))
    }

    if (generic == "gamma") {
        return(apply_delayed_gamma(seed))
    }

    if (generic == "lgamma") {
        return(apply_delayed_lgamma(seed))
    }

    if (generic == "round") {
        if (envir$digits != 0) {
            return("only 'digits = 0' are supported for delayed 'round'")
        }
        return(apply_delayed_round(seed))
    }

    log.base.support <- c(log2=2, log10=10)
    if (generic %in% names(log.base.support)) {
        return(apply_delayed_log(seed, log.base.support[[generic]]))
    }

    NULL
}

.unary_Ops <- function(seed, OP) {
    envir <- environment(OP)
    generic <- envir$`.Generic`

    if (is.null(generic)) {
        return(NULL)
    }

    if (generic == "!") {
        return(apply_delayed_boolean_not(seed))
    }

    if (!(generic %in% supported.Ops)) {
        return(NULL)
    }

    e1 <- envir$e1
    e2 <- envir$e2

    if (missing(e2)) {
        if (generic == "+") {
            return(seed)
        } else if (generic == "-") {
            return(apply_delayed_associative_arithmetic(seed, -1, TRUE, "*"))
        } else {
            stop("second argument can only be missing for unary '+' or '-'")
        }
    } 

    # i.e., is the operation applied to the right of the seed?
    right <- is(e1, "DelayedArray") 
    val <- if (right) e2 else e1

    # Just pretending that we're applying by rows; this should never happen,
    # as vector operations are handled by DelayedUnaryIsoOpWithArgs.
    output <- .apply_delayed_unary_ops(seed, generic, val, right, row=TRUE)
    stopifnot(!is.null(output))
    output
}

#' @export
setMethod("initializeCpp", "DelayedUnaryIsoOpStack", function(x, ...) {
    seed <- initializeCpp(x@seed, ...)

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

        if (!status || is.character(seed)) {
            msg <- if (is.character(seed)) seed else paste0("unsupported function in '<", class(x), ">@OPS[[", i, "]]'")
            warning(msg, ", falling back to an unknown matrix")
            return(initialize_unknown_matrix(x))
        }
    }

    seed
})

####################################################################################
####################################################################################

#' @export
setMethod("initializeCpp", "DelayedNaryIsoOp", function(x, ...) {
    if (length(x@seeds) != 2) {
        stop("expected exactly two seeds for 'DelayedNaryIsoOp'")
    }
    if (length(x@Rargs)) {
        stop("expected no additional right arguments for 'DelayedNaryIsoOp'")
    }

    left <- initializeCpp(x@seeds[[1]], ...)
    right <- initializeCpp(x@seeds[[2]], ...)

    # Figuring out the identity of the operation.
    chosen <- NULL
    for (p in supported.Ops) {
        if (identical(x@OP, get(p, envir=baseenv()))) {
            chosen <- p
            break
        }
    }

    if (is.null(chosen)) {
        warning("unknown operation in '<", class(x)[1], ">@OP', falling back to an unknown matrix")
        return(initialize_unknown_matrix(x))
    }

    output <- apply_delayed_binary_operation(left, right, chosen)
    stopifnot(!is.null(output))
    output
})
