use_beachmat <- function() {
    if (!requireNamespace("usethis", quietly = TRUE)) {
        stop("usethis package required for this function.")
    }
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
    usethis::use_directory("src")
    usethis::use_git_ignore(c("*.o", "*.so", "*.dll"), "src")
    if (usethis:::uses_roxygen()) {
        usethis:::todo(
            "Include the following roxygen tags somewhere in your package")
        usethis:::code_block(
            paste0(
                "#' @useDynLib ",
                usethis:::project_name(),
                ", .registration = TRUE"),
            "#' @importFrom Rcpp sourceCpp",
            "NULL",
            copy = TRUE)
    }
    else {
        usethis:::todo("Include the following directives in your NAMESPACE")
        usethis:::code_block(
            paste0(
                "useDynLib('",
                usethis:::project_name(),
                "', .registration = TRUE)"),
            "importFrom('Rcpp', 'sourceCpp')",
            copy = TRUE)
        usethis:::edit_file(usethis:::proj_path("NAMESPACE"))
    }
    usethis:::write_over(
        base_path = usethis:::scoped_path("project", "src"),
        path = "Makevars",
        contents = paste0(
            "BEACHMAT_LIBS=`echo 'beachmat::pkgconfig(\"PKG_LIBS\")'|\\\n",
            "      \"${R_HOME}/bin/R\" --vanilla --slave`\n",
            "PKG_LIBS=$(BEACHMAT_LIBS)\n"))
    usethis:::write_over(
        base_path = usethis:::scoped_path("project", "src"),
        path = "Makevars.win",
        contents = paste0(
            "BEACHMAT_LIBS=$(shell echo 'beachmat::pkgconfig(\"PKG_LIBS\")'|\\\n",
            "      \"${R_HOME}/bin/R\" --vanilla --slave)\n",
            "PKG_LIBS=$(BEACHMAT_LIBS)\n"))
    usethis:::todo(
        "If your package needs to add to the $PKG_LIBS variable ",
        "in Makevars and Makevars.win, please consult the ",
        "'Linking' vignette available via:")
    usethis:::code_block(
        "vignette(\"linking\", package = \"beachmat\")",
        copy = FALSE)
    usethis:::todo("Run document()")
}
