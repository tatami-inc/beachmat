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
        usethis:::todo("Include the following roxygen tags somewhere in your ",
                       "package")
        usethis:::code_block(paste0("#' @useDynLib ",
                                    usethis:::project_name(),
                                    ", .registration = TRUE"),
                             "#' @importFrom Rcpp sourceCpp", "NULL",
                             copy = FALSE)
    }
    else {
        usethis:::todo("Include the following directives in your NAMESPACE")
        usethis:::code_block(paste0("useDynLib('",
                                    usethis:::project_name(),
                                    "', .registration = TRUE)"),
                             "importFrom('Rcpp', 'sourceCpp')",
                             copy = FALSE)
        usethis:::edit_file(usethis:::proj_path("NAMESPACE"))
    }
    usethis:::todo("Add the following directives in your src/Makevars")
    usethis:::code_block(paste0("BEACHMAT_LIBS=`echo ",
                                "'beachmat::pkgconfig(\"PKG_LIBS\")'|\\\n",
                                "      \"${R_HOME}/bin/R\" --vanilla --slave`"),
                         "PKG_LIBS=$(BEACHMAT_LIBS)",
                         copy = FALSE)
    usethis:::edit_file(usethis:::scoped_path("project", "src", "Makevars"))
    usethis:::todo("Add the following directives in your src/Makevars.win ",
                   "(be careful of subtle differences to src/Makevars)")
    usethis:::code_block(paste0("BEACHMAT_LIBS=$(shell echo ",
                                "'beachmat::pkgconfig(\"PKG_LIBS\")'|\\\n",
                                "      \"${R_HOME}/bin/R\" --vanilla --slave)"),
                         "PKG_LIBS=$(BEACHMAT_LIBS)",
                         copy = FALSE)
    usethis:::edit_file(usethis:::scoped_path("project", "src", "Makevars.win"))
    usethis:::todo("If your package needs to add to the $PKG_LIBS variable, ",
                   "please consult the 'Linking' vignette available via:")
    usethis:::code_block("vignette(\"linking\", package = \"beachmat\")",
                         copy = FALSE)
    usethis:::todo("Restart R for changes to take effect")
    usethis:::todo("Run document()")
}
