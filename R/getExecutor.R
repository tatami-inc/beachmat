#' Get the parallel executor
#'
#' Get the executor object for safe execution of R code in parallel sections.
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
