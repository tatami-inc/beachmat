#' Check the in-memory cache for matrix instances
#'
#' Check the in-memory cache for a pre-existing initialized C++ object, and initialize it if it does not exist.
#' This is typically used in \code{\link{initializeCpp}} methods of file-backed representations to avoid redundant reads of the entire matrix.
#'
#' @param namespace String containing the namespace, typically the name of the package implementing the method.
#' @param key String containing the key for a specific matrix instance.
#' @param fun Function that accepts no arguments and returns an external pointer like those returned by \code{\link{initializeCpp}}.
#'
#' @return For \code{checkMemoryCache}, the output of \code{fun} (possibly from an existing cache) is returned.
#'
#' For \code{flushMemoryCache}, all existing cached objects are removed and \code{NULL} is invisibly returned.
#'
#' @details
#' For representations where data extraction is costly (e.g., from file), \code{\link{initializeCpp}} methods may provide a \code{memorize=} option.
#' Setting this to \code{TRUE} will load the entire matrix into memory, effectively paying a one-time up-front cost to improve efficiency for downstream operations that pass through the matrix multiple times.
#'
#' If this option is provided, \code{initializeCpp} methods are expected to cache the in-memory instance using \code{checkMemoryCache}.
#' This ensures that all subsequent calls to the same \code{initializeCpp} method will return the same instance, avoiding redundant memory loads when the same matrix is used in multiple functions.
#'
#' Of course, this process saves time at the expense of increased memory usage.
#' If too many instances are being cached, they can be cleared from memory using the \code{flushMemoryCache} function.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking up a class with some kind of uniquely identifying aspect.
#' setClass("UnknownMatrix", slots=c(contents="dgCMatrix", uuid="character"))
#' X <- new("UnknownMatrix", 
#'     contents=Matrix::rsparsematrix(10, 10, 0.1), 
#'     uuid=as.character(sample(1e8, 1)))
#'
#' # Defining our initialization method.
#' setMethod("initializeCpp", "UnknownMatrix", function(x, ..., memorize=FALSE) {
#'     if (memorize) {
#'         checkMemoryCache("my_package", x@uuid, function() initializeCpp(x@contents))
#'     } else {
#'         initializeCpp(x@contents)
#'     }
#' })
#'
#' # Same pointer is returned multiple times.
#' initializeCpp(X, memorize=TRUE)
#' initializeCpp(X, memorize=TRUE)
#' 
#' # Flushing the cache.
#' flushMemoryCache()
#'
#' @name checkMemoryCache
NULL

memory.cache <- new.env()
memory.cache$contents <- list()

#' @export
#' @rdname checkMemoryCache
flushMemoryCache <- function() {
    memory.cache$contents <- list()
    gc()
    invisible(NULL)
}

#' @export
#' @rdname checkMemoryCache
checkMemoryCache <- function(namespace, key, fun) {
    ns <- memory.cache$contents[[namespace]]
    if (is.null(ns)) {
        ns <- list()
    }

    if (key %in% names(ns)) {
        return(ns[[key]])
    }

    ns[[key]] <- fun()
    memory.cache$contents[[namespace]] <- ns
    return(ns[[key]])
}
