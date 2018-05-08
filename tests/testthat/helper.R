proj <- new.env(parent = emptyenv())

# Copied from
# https://github.com/r-lib/usethis/blob/master/tests/testthat/helper.R

## putting `pattern` in the package or project name is part of our strategy for
## suspending the nested project check during testing
pattern <- "aaa"

scoped_temporary_thing <- function(dir = tempfile(pattern = pattern),
                                   env = parent.frame(),
                                   rstudio = FALSE,
                                   thing = c("package", "project")) {
    if (requireNamespace("withr", quietly = TRUE) &&
        requireNamespace("usethis", quietly = TRUE)) {
        thing <- match.arg(thing)
        old <- proj$cur
        # Can't schedule a deferred project reset if calling this from the R
        # console, which is useful when developing tests
        if (identical(env, globalenv())) {
            usethis:::todo(
                "Switching to a temporary project! To restore current project:\n",
                "proj_set(\"", old, "\")"
            )
        } else {
            withr::defer(proj_set(old), envir = env)
        }

        switch(
            thing,
            package = capture_output(
                usethis::create_package(dir, rstudio = rstudio, open = FALSE)),
            project = capture_output(
                usethis::create_project(dir, rstudio = rstudio, open = FALSE))
        )
        invisible(dir)
    }
}

scoped_temporary_project <- function(dir = tempfile(pattern = pattern),
                                     env = parent.frame(),
                                     rstudio = FALSE) {
    scoped_temporary_thing(dir, env, rstudio, "project")
}

scoped_temporary_package <- function(dir = tempfile(pattern = pattern),
                                     env = parent.frame(),
                                     rstudio = FALSE) {
    scoped_temporary_thing(dir, env, rstudio, "package")
}
