#' Use beachmat
#'
#' Creates \code{src/} and adds needed packages to \code{DESCRIPTION}.
#'
#' For more details see the 'Linking' vignette, available by running
#' \code{vignette("linking", package = "beachmat")}.
#' @importFrom usethis use_directory use_git_ignore
#' @export
use_beachmat <- function() {
  usethis:::check_is_package("use_beachmat()")
  usethis:::use_dependency("Rcpp", "LinkingTo")
  usethis:::use_dependency("Rhdf5lib", "LinkingTo")
  usethis:::use_dependency("beachmat", "LinkingTo")
  usethis:::use_dependency("Rcpp", "Imports")
  usethis:::use_dependency("beachmat", "Suggests")
  usethis:::use_dependency("HDF5Array", "Suggests")
  usethis:::use_dependency("DelayedArray", "Suggests")
  usethis:::use_dependency("Matrix", "Suggests")
  usethis:::use_description_field("SystemRequirements", "C++11")
  use_directory("src")
  use_git_ignore(c("*.o", "*.so", "*.dll"), "src")
  if (usethis:::uses_roxygen()) {
    usethis:::todo(
      "Include the following roxygen tags somewhere in your package")
    usethis:::code_block(paste0(
      "#' @useDynLib ", usethis:::project_name(), ", .registration = TRUE"),
      "#' @importFrom Rcpp sourceCpp", "NULL")
  }
  else {
    usethis:::todo("Include the following directives in your NAMESPACE")
    usethis:::code_block(paste0(
      "useDynLib('", usethis:::project_name(), "', .registration = TRUE)"),
      "importFrom('Rcpp', 'sourceCpp')")
    usethis:::edit_file(usethis:::proj_path("NAMESPACE"))
  }
  usethis:::todo("Run document()")
}
