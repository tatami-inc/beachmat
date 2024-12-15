#' Tatami utilities
#'
#' Utility functions that directly operate on the pointers produced by \code{\link{initializeCpp}}.
#' Some of these are used internally by \code{initializeCpp} methods operating on \pkg{DelayedArray} classes.
#'
#' @param x A pointer produced by \code{\link{initializeCpp}}.
#' @param xs A list of pointers produced by \code{\link{initializeCpp}}.
#' All matrices should have the same number of rows (if \code{by.row=FALSE}) or columns (otherwise).
#' @param by.row Logical scalar indicating whether to apply the operation on the rows.
#' \itemize{
#' \item For \code{tatami.bind}, this will combine the matrices by rows,
#' i.e., the output matrix has a number of rows equal to the sum of the number of rows in \code{xs}.
#' \item For \code{tatami.subset}, this will subset the matrix by row.
#' \item For \code{tatami.arith}, \code{tatami.compare} and \code{tatami.logic} with a vector \code{val},
#' the vector should have length equal to the number of rows.k
#' }
#' @param op String specifying the operation to perform.
#' \itemize{
#' \item For \code{tatami.arith}, this should be one of the operations in \link{Arith}.
#' \item For \code{tatami.compare}, this should be one of the operations in \link{Compare}.
#' \item For \code{tatami.logic}, this should be one of the operations in \link{Logic}.
#' \item For \code{tatami.math}, this should be one of the operations in \link{Math}.
#' \item For \code{tatami.binary}, this may be any operation in \link{Arith}, \link{Compare} or \link{Logic}.
#' }
#' @param val For \code{tatami.arith}, \code{tatami.compare} and \code{tatami.logic}, the value to be used in the operation specified by \code{op}. 
#' This may be a:
#' \itemize{
#' \item Numeric scalar, which is used in the operation for all entries of the matrix.
#' \item Numeric vector of length equal to the number of rows, where each value is used in the operation with the corresponding row when \code{by.row=TRUE}.
#' \item Numeric vector of length equal to the number of column, where each value is used with the corresponding column when \code{by.row=FALSE}.
#' }
#'
#' For \code{tatami.multiply}, the value to be used in the matrix multiplication.
#' This may be a:
#' \itemize{
#' \item Numeric vector of length equal to the number of columns of \code{x} (if \code{right=FALSE}) or rows (otherwise).
#' \item Numeric matrix with number of rows equal to the number of columns of \code{x} (if \code{right=FALSE}) or rows (otherwise).
#' \item Pointer produced by \code{\link{initializeCpp}}, 
#' referencing a matrix with number of rows equal to the number of columns of \code{x} (if \code{right=FALSE}) or rows (otherwise).
#' }
#' @param right For \code{tatami.arith} and \code{tatami.compare}, 
#' a logical scalar indicating that \code{val} is on the right-hand side of the operation.
#'
#' For \code{tatami.multiply}, a logical scalar indicating that \code{val} is on the right-hand side of the multiplication.
#' @param subset Integer vector containing the subset of interest.
#' These should be 1-based row or column indices depending on \code{by.row}.
#' @param y A pointer produced by \code{\link{initializeCpp}},
#' referencing a matrix of the same dimensions as \code{x}.
#' @param base Numeric scalar specifying the base of the log-transformation.
#' @param i Integer scalar containing the 1-based index of the row (for \code{tatami.row}) or column (for \code{tatami.column}) of interest. 
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return 
#' For \code{tatami.dim}, an integer vector containing the dimensions of the matrix.
#'
#' For \code{tatami.is.sparse}, a logical scalar indicating whether the matrix is sparse.
#'
#' For \code{tatami.prefer.rows}, a logical scalar indicating whether the matrix prefers iteration by row.
#'
#' For \code{tatami.row} or \code{tatami.column}, a numeric vector containing the contents of row or column \code{i}, respectively.
#'
#' For \code{tatami.row.sums} or \code{tatami.column.sums}, a numeric vector containing the row or column sums, respectively.
#'
#' For \code{tatami.row.nan.counts} or \code{tatami.column.nan.counts}, a numeric vector containing the number of NaNs in each row or column, respectively.
#'
#' For \code{tatami.realize}, a numeric matrix or \linkS4class{dgCMatrix} with the matrix contents.
#' The exact class depends on whether \code{x} refers to a sparse matrix. 
#' 
#' For \code{tatami.multiply}, a numeric matrix containing the matrix product of \code{x} and \code{other}.
#' 
#' For all other functions, a new pointer to a matrix with the requested operations applied to \code{x} or \code{xs}.
#'
#' @author Aaron Lun
#' @examples
#' x <- Matrix::rsparsematrix(1000, 100, 0.1)
#' ptr <- initializeCpp(x)
#' tatami.dim(ptr)
#' tatami.row(ptr, 1)
#'
#' rounded <- tatami.round(ptr)
#' tatami.row(rounded, 1)
#'
#' @name tatami-utils
NULL

#' @export
#' @rdname tatami-utils
tatami.bind <- function(xs, by.row) { 
    apply_delayed_bind(xs, by.row)
}

#' @export
#' @rdname tatami-utils
tatami.transpose <- function(x) {
    apply_delayed_transpose(x)
}

#' @export
#' @rdname tatami-utils
tatami.subset <- function(x, subset, by.row) {
    apply_delayed_subset(x, subset, by.row)
}

#' @export
#' @rdname tatami-utils
tatami.arith <- function(x, op, val, by.row, right) {
    if (op %in% supported.Arith1) {
        apply_delayed_associative_arithmetic(x, val, by.row, op)
    } else if (op %in% supported.Arith2) {
        apply_delayed_nonassociative_arithmetic(x, val, right, by.row, op)
    } else {
        stop("unknown operation '", op, "'")
    }
}

#' @export
#' @rdname tatami-utils
tatami.compare <- function(x, op, val, by.row, right) {
    op <- match.arg(op, supported.Compare)
    if (!right) { # need to flip the operation if the argument is not on the right.
        op <- reverse.Compare[[op]]
    }
    apply_delayed_comparison(x, val, by.row, op)
}

#' @export
#' @rdname tatami-utils
tatami.logic <- function(x, op, val, by.row) {
    op <- match.arg(op, supported.Logic)
    apply_delayed_boolean(x, val, by.row, op)
}

#' @export
#' @rdname tatami-utils
tatami.round <- function(x) {
    apply_delayed_round(x)
}

#' @export
#' @rdname tatami-utils
tatami.log <- function(x, base) {
    apply_delayed_log(x, base)
}

#' @export
#' @rdname tatami-utils
tatami.math <- function(x, op) {
    apply_delayed_unary_math(x, op)
}

#' @export
#' @rdname tatami-utils
tatami.not <- function(x) {
    apply_delayed_boolean_not(x)
}

#' @export
#' @rdname tatami-utils
tatami.binary <- function(x, y, op) {
    op <- match.arg(op, supported.Ops)
    apply_delayed_binary_operation(x, y, op)
}

#' @export
#' @rdname tatami-utils
tatami.dim <- function(x) {
    tatami_dim(x)
}

#' @export
#' @rdname tatami-utils
tatami.row <- function(x, i) {
    tatami_row(x, i)
}

#' @export
#' @rdname tatami-utils
tatami.column <- function(x, i) {
    tatami_column(x, i)
}

#' @export
#' @rdname tatami-utils
tatami.row.sums <- function(x, num.threads) {
    tatami_row_sums(x, num.threads)
}

#' @export
#' @rdname tatami-utils
tatami.column.sums <- function(x, num.threads) {
    tatami_column_sums(x, num.threads)
}

#' @export
#' @rdname tatami-utils
tatami.row.nan.counts <- function(x, num.threads) {
    tatami_row_nan_counts(x, num.threads)
}

#' @export
#' @rdname tatami-utils
tatami.column.nan.counts <- function(x, num.threads) {
    tatami_column_nan_counts(x, num.threads)
}

#' @export
#' @rdname tatami-utils
tatami.is.sparse <- function(x) {
    tatami_is_sparse(x)
}

#' @export
#' @rdname tatami-utils
tatami.prefer.rows <- function(x) {
    tatami_prefer_rows(x)
}

#' @export
#' @rdname tatami-utils
tatami.realize <- function(x, num.threads) {
    tatami_realize(x, num.threads)
}

#' @export
#' @rdname tatami-utils
tatami.multiply <- function(x, val, right, num.threads) {
    if (is.atomic(val)) {
        if (is.null(dim(val))) {
            tatami_multiply_vector(x, val, right=right, num_threads=num.threads)
        } else if (!right) {
            t(tatami_multiply_columns(x, t(val), right=right, num_threads=num.threads))
        } else {
            tatami_multiply_columns(x, val, right=right, num_threads=num.threads)
        }
    } else {
        tatami_multiply_matrix(x, val, right=right, num_threads=num.threads)
    }
}
