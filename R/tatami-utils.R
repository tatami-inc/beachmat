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
#' @param row For \code{tatami.get}, a boolean indicating whether to extract the \code{i}-th row.
#' If \code{FALSE}, the \code{i}-th column is extracted instead.
#'
#' For \code{tatami.sums}, \code{tatami.sums.by.group}, \code{tatami.medians}, etc., a boolean indicating whether to compute row-wise statistics.
#' If \code{FALSE}, column-wise statistics are computed instead.
#' @param group Integer vector of length equal to the number of columns (if \code{row = TRUE}) or rows (otherwise),
#' containing the group assignment for each column and row, respectively.
#' Assignments should lie in \code{[1, N]} where \code{N} is the total number of groups.
#' @param subset Integer vector containing the subset of interest.
#' These should be 1-based row or column indices depending on \code{by.row}.
#' @param y A pointer produced by \code{\link{initializeCpp}},
#' referencing a matrix of the same dimensions as \code{x}.
#' @param base Numeric scalar specifying the base of the log-transformation.
#' @param i Integer scalar containing the 1-based index of the row (for \code{row=TRUE}) or column (otherwise) of interest. 
#' This should be in \code{[1, D]} where \code{D} is the total number of rows or columns, respectively, in \code{x}.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return 
#' For \code{tatami.dim}, an integer vector containing the dimensions of the matrix.
#'
#' For \code{tatami.is.sparse}, a logical scalar indicating whether the matrix is sparse.
#'
#' For \code{tatami.prefer.rows}, a logical scalar indicating whether the matrix prefers iteration by row.
#'
#' For \code{tatami.get}, a numeric vector containing the contents of row or column \code{i}, respectively.
#'
#' For \code{tatami.realize}, a numeric matrix or dgCMatrix with the matrix contents.
#' The exact class depends on whether \code{x} refers to a sparse matrix. 
#' 
#' For \code{tatami.multiply}, a numeric matrix containing the matrix product of \code{x} and \code{other}.
#'
#' For \code{tatami.sums}, a numeric vector containing the row or column sums, respectively.
#'
#' For \code{tatami.sums.by.group}, a numeric matrix is returned.
#' \itemize{
#' \item If \code{row=TRUE}, its number of rows is equal to that of \code{x} and its number of columns is equal to the number of groups.
#' Each column corresponds to a group and contains the row sums across all columns of \code{x} assigned to that group.
#' \item If \code{row=FALSE}, its number of columns is equal to that of \code{x} and its number of rows is equal to the number of groups.
#' Each row corresponds to a group and contains the column sums across all columns of \code{x} assigned to that group.
#' }
#'
#' For \code{tatami.medians}, a numeric vector containing the row or column medians, respectively.
#'
#' For \code{tatami.nan.counts}, a numeric vector containing the number of NaNs in each row or column, respectively.
#' 
#' For all other functions, a new pointer to a matrix with the requested operations applied to \code{x} or \code{xs}.
#'
#' @aliases
#' tatami.row.medians
#' tatami.column.medians
#' tatami.row.sums
#' tatami.column.sums
#' tatami.row.nan.counts
#' tatami.column.nan.counts
#'
#' @author Aaron Lun
#' @examples
#' x <- Matrix::rsparsematrix(1000, 100, 0.1)
#' ptr <- initializeCpp(x)
#' tatami.dim(ptr)
#' tatami.get(ptr, 1, row=TRUE)
#'
#' rounded <- tatami.round(ptr)
#' tatami.get(rounded, 1, row=TRUE)
#'
#' tatami.sums(ptr, row=FALSE, num.threads=2)
#' tatami.medians(ptr, row=FALSE, num.threads=2)
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
tatami.get <- function(x, i, row) {
    tatami_get(x, i, row)
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
            tatami_multiply_vector(x, val, right=right, threads=num.threads)
        } else if (!right) {
            t(tatami_multiply_columns(x, t(val), right=right, threads=num.threads))
        } else {
            tatami_multiply_columns(x, val, right=right, threads=num.threads)
        }
    } else {
        tatami_multiply_matrix(x, val, right=right, threads=num.threads)
    }
}

#' @export
#' @rdname tatami-utils
tatami.sums <- function(x, row, num.threads) {
    tatami_sums(x, row, num.threads)
}

#' @export
#' @rdname tatami-utils
tatami.sums.by.group <- function(x, group, row, num.threads) {
    tatami_sums_by_group(x, group, row, num.threads) 
}

#' @export
#' @rdname tatami-utils
tatami.medians <- function(x, row, num.threads) {
    tatami_medians(x, row, num.threads)
}

#' @export
#' @rdname tatami-utils
tatami.nan.counts <- function(x, row, num.threads) {
    tatami_nan_counts(x, row, num.threads)
}

##### For compatibility only. #######

#' @export
tatami.row <- function(x, i) {
    tatami.get(x, i, row=TRUE)
}

#' @export
tatami.column <- function(x, i) {
    tatami.get(x, i, row=FALSE)
}

#' @export
tatami.row.sums <- function(x, num.threads) {
    tatami.sums(x, row=TRUE, num.threads=num.threads)
}

#' @export
tatami.column.sums <- function(x, num.threads) {
    tatami.sums(x, row=FALSE, num.threads=num.threads)
}

#' @export
tatami.row.medians <- function(x, num.threads) {
    tatami.medians(x, row=TRUE, num.threads=num.threads)
}

#' @export
tatami.column.medians <- function(x, num.threads) {
    tatami.medians(x, row=FALSE, num.threads=num.threads)
}

#' @export
tatami.row.nan.counts <- function(x, num.threads) {
    tatami.nan.counts(x, row=TRUE, num.threads=num.threads)
}

#' @export
tatami.column.nan.counts <- function(x, num.threads) {
    tatami.nan.counts(x, row=FALSE, num.threads=num.threads)
}
