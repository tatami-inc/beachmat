#' Get the parallel executor
#'
#' Get the executor object for safe execution of R code in parallel sections.
#' This should be set by \code{Rtatami::set_executor()} in the \code{.onLoad} function of downstream packages.
#'
#' @return An external pointer to be passed to \code{Rtatami::set_executor}.
#'
#' @author Aaron Lun
#' @examples
#' getExecutor()
#' 
#' @export
getExecutor <- function() {
    get_executor()
}
