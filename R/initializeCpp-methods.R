is_class_package <- function(x, package, classes) {
    if (isNamespaceLoaded(package)) {
        for (cls in classes) {
            attr(cls, "package") <- package
            if (is(x, cls)) {
                return(TRUE)
            }
        }
    }

    FALSE
}

#' @export
setMethod("initializeCpp", "ANY", function(x, ...) {
    if (is_class_package(x, "HDF5Array", c("HDF5ArraySeed", "H5SparseMatrixSeed"))) {
        # Automatically use beachmat.hdf5 if it's available.  We check that
        # beachmat.hdf5 is not already loaded to avoid infinite loops in case
        # beachmat.hdf5 does NOT support the listed classes (so requiring the
        # namespace would not add any more methods for meaningful dispatch).
        if (!isNamespaceLoaded("beachmat.hdf5") && requireNamespace("beachmat.hdf5", quietly=TRUE)) {
            return(initializeCpp(x, ...))
        }
    }

    if (is_class_package(x, "TileDBArray", "TileDBArraySeed")) {
        # Same for the TileDB matrices.
        if (!isNamespaceLoaded("beachmat.tiledb") && requireNamespace("beachmat.tiledb", quietly=TRUE)) {
            return(initializeCpp(x, ...))
        }
    }

    if (is_class_package(x, "alabaster.matrix", c("WrapperArraySeed", "ReloadedArraySeed"))) {
        # Pass-through some known no-op matrices from alabaster.matrix.
        return(initializeCpp(x@seed, ...))
    }

    # We also provide a lightweight hook to support skipping of classes in
    # other packages (mostly legacy GNE-internal gunk) that don't want to
    # import beachmat just to define a no-op initializeCpp() method.
    for (cls in getOption("beachmat.noop.wrappers", character(0))) {
        if (is(x, cls)) {
            return(initializeCpp(x@seed, ...))
        }
    }

    initialize_unknown_matrix(x)
})

#' @export
setMethod("initializeCpp", "externalptr", function(x, ...) x)

#' @export
setMethod("initializeCpp", "matrix", function(x, .check.na = TRUE, ...) initialize_dense_matrix(x, nrow(x), ncol(x), check_na=.check.na))

#' @export
setMethod("initializeCpp", "dgeMatrix", function(x, ...) initialize_dense_matrix_from_vector(x@x, nrow(x), ncol(x), check_na=FALSE))

#' @export
setMethod("initializeCpp", "lgeMatrix", function(x, .check.na = TRUE, ...) initialize_dense_matrix_from_vector(x@x, nrow(x), ncol(x), check_na=.check.na))

####################################################################################
####################################################################################

#' @export
#' @import Matrix
setMethod("initializeCpp", "dgCMatrix", function(x, ...) initialize_sparse_matrix(x@x, x@i, x@p, nrow(x), ncol(x), byrow=FALSE, check_na=FALSE))

#' @export
setMethod("initializeCpp", "dgRMatrix", function(x, ...) initialize_sparse_matrix(x@x, x@j, x@p, nrow(x), ncol(x), byrow=TRUE, check_na=FALSE))

#' @export
#' @import Matrix
setMethod("initializeCpp", "lgCMatrix", function(x, .check.na = TRUE, ...) initialize_sparse_matrix(x@x, x@i, x@p, nrow(x), ncol(x), byrow=FALSE, check_na=.check.na))

#' @export
setMethod("initializeCpp", "lgRMatrix", function(x, .check.na = TRUE, ...) initialize_sparse_matrix(x@x, x@j, x@p, nrow(x), ncol(x), byrow=TRUE, check_na=.check.na))

####################################################################################
####################################################################################

#' @export
#' @import DelayedArray 
setMethod("initializeCpp", "DelayedMatrix", function(x, ...) {
    tryCatch(
        initializeCpp(x@seed, ...),
        error=function(e) {
            warning(gsub("\\s+$", " ", e$message), ", falling back to an unknown matrix")
            initialize_unknown_matrix(x)
        }
    )
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
#' @importClassesFrom DelayedArray ConstantMatrix
setMethod("initializeCpp", "ConstantArraySeed", function(x, ...) {
    initialize_constant_matrix(x@dim[1], x@dim[2], x@value)
})

####################################################################################
####################################################################################

#' @export
#' @importClassesFrom SparseArray SVT_SparseMatrix
setMethod("initializeCpp", "SVT_SparseMatrix", function(x, .check.na = TRUE, ...) {
    initialize_SVT_SparseMatrix(nr=nrow(x), nc=ncol(x), x, check_na = .check.na)
})

####################################################################################
####################################################################################

supported.Arith1 <- c("+", "*")
supported.Arith2 <- c("-", "/", "^", "%/%", "%%")
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
        stop("unknown operation in '<", class(x)[1], ">@OP'")
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

    if (generic == "round") {
        if (envir$digits != 0) {
            stop("only 'digits = 0' are supported for delayed 'round'")
        }
        return(apply_delayed_round(seed))
    }

    log.base.support <- c(log2=2, log10=10)
    if (generic %in% names(log.base.support)) {
        return(apply_delayed_log(seed, log.base.support[[generic]]))
    }

    apply_delayed_unary_math(seed, generic)
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

    if (generic == "type<-") {
        target.type <- envir$value
        src.type <- type(envir$x)

        ranking <- c("logical", "raw", "integer", "double")
        target.i <- which(ranking == target.type)
        src.i <- which(ranking == src.type)
        if (length(target.i) == 0L || length(src.i) == 0L) {
            return(NULL)
        }

        if (target.i >= src.i) {
            return(seed)
        }

        if (target.type == "integer") {
            # source type must be a wider type, so we truncate.
            return(apply_delayed_unary_math(seed, "trunc"))
        } else if (target.type == "logical") {
            # TODO: add a better 'not-not' method to coerce values to binary.
            return(apply_delayed_boolean(seed, val=TRUE, row=TRUE, op="&")) 
        } else {
            # Raw not supported yet.
            return(NULL)
        }
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

        if (!status) {
            stop("unsupported function in '<", class(x), ">@OPS[[", i, "]]'")
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
        stop("unknown operation in '<", class(x)[1], ">@OP'")
    }

    output <- apply_delayed_binary_operation(left, right, chosen)
    stopifnot(!is.null(output))
    output
})
