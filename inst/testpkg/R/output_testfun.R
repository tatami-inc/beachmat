# Creating functions to check each set of matrix inputs.

###############################

.check_output_mat <- function(FUN, ..., class.out, cxxfun) { 
    for (i in 1:3) { 
        test.mat <- FUN(...)

        # Specifying order.
        if (i==1L) {
            ranges <- list(forward=seq_len(ncol(test.mat))-1L,
                           random=sample(ncol(test.mat))-1L)
            ranges$reverse <- rev(ranges$forward)
        } else if (i==2L) {
            ranges <- list(forward=seq_len(nrow(test.mat))-1L,
                           random=sample(nrow(test.mat))-1L)
            ranges$reverse <- rev(ranges$forward)
        } else {
            ranges <- list(integer(0)) # doesn't matter.
        }

        # We should get the same results, regardless of the order.
        for (ordering in ranges) { 
            out <- .Call(cxxfun, test.mat, i, ordering)
    
            if (is.matrix(test.mat)) { 
                testthat::expect_true(is.matrix(out[[1]]))
            } else {
                class.out <- class(test.mat)
                testthat::expect_s4_class(out[[1]], class.out)

                if (class.out=="HDF5Matrix") { 
                    testthat::expect_equal(out[[1]]@seed@first_val, as.vector(out[[1]][1,1]))
                    testthat::expect_identical(dim(out[[1]]), dim(test.mat))
                    testthat::expect_true(path(out[[1]])!=path(test.mat))
                }

                out[[1]] <- as.matrix(out[[1]])
            }      

            # Checking that the set-and-get values are the same.
            ref <- as.matrix(test.mat)
            ref2 <- matrix(out[[2]], nrow(ref), ncol(ref), byrow=(i==2L))
            testthat::expect_identical(ref, ref2)

            # Reordering 'ref' to match the expected output. 
            if (i==1L) {
                ref[,ordering+1L] <- ref
            } else if (i==2L) {
                ref[ordering+1L,] <- ref
            }
            testthat::expect_identical(ref, out[[1]])
        }
    }
    return(invisible(NULL))
}
   
check_integer_output_mat <- function(FUN, ...) {
    .check_output_mat(FUN=FUN, ..., cxxfun=cxx_test_integer_output)
} 

check_numeric_output_mat <- function(FUN, ...) {
    .check_output_mat(FUN=FUN, ..., cxxfun=cxx_test_numeric_output)
} 

check_logical_output_mat <- function(FUN, ...) {
    .check_output_mat(FUN=FUN, ..., cxxfun=cxx_test_logical_output)
} 

check_character_output_mat <- function(FUN, ...) {
    .check_output_mat(FUN=FUN, ..., cxxfun=cxx_test_character_output)
} 

###############################

.check_output_slice <- function(FUN, ..., by.row, by.col, cxxfun, fill) { 
    rx <- range(by.row)
    ry <- range(by.col)

    for (it in 1:2) {
        test.mat <- FUN(...)      
        out <- .Call(cxxfun, test.mat, it, rx, ry)

        if (is.matrix(test.mat)) { 
            testthat::expect_true(is.matrix(out[[1]]))
        } else {
            testthat::expect_s4_class(out[[1]], class(test.mat))
            out[[1]] <- as.matrix(out[[1]])
        }

        # Checking the standard fill.
        test <- out[[1]]
        ref <- as.matrix(test.mat)
        testthat::expect_identical(ref[by.row, by.col], test[by.row, by.col]) 
        test[by.row,by.col] <- fill
        testthat::expect_true(all(test==fill))

        # Checking the set-and-get.
        ref2 <- matrix(out[[2]], diff(rx)+1L, diff(ry)+1L, byrow=(it==2))
        testthat::expect_identical(ref[by.row, by.col], ref2) 
    }
    return(invisible(NULL))
}


check_integer_output_slice <- function(FUN, ..., by.row, by.col) {
    .check_output_slice(FUN=FUN, ..., by.row=by.row, by.col=by.col, 
                        cxxfun=cxx_test_integer_output_slice, fill=0L)
} 

check_numeric_output_slice <- function(FUN, ..., by.row, by.col) {
    .check_output_slice(FUN=FUN, ..., by.row=by.row, by.col=by.col, 
                        cxxfun=cxx_test_numeric_output_slice, fill=0)
} 

check_logical_output_slice <- function(FUN, ..., by.row, by.col) {
    .check_output_slice(FUN=FUN, ..., by.row=by.row, by.col=by.col, 
                        cxxfun=cxx_test_logical_output_slice, fill=FALSE)
} 

check_character_output_slice <- function(FUN, ..., by.row, by.col) {
    .check_output_slice(FUN=FUN, ..., by.row=by.row, by.col=by.col, 
                        cxxfun=cxx_test_character_output_slice, fill="")
}

###############################

.check_output_indexed <- function(FUN, ..., N, cxxfun, fill) { 
    for (n in N) {
        for (it in 1:1) {
            test.mat <- FUN(...)
            all.values <- as.matrix(test.mat)
            ref <- matrix(fill, nrow(test.mat), ncol(test.mat))

            if (it==1L) {
                # Constructing a list of row index vectors for a sampled set of columns.
                c <- sample(ncol(test.mat), n, replace=TRUE)
                subr <- vector("list", n)
                for (i in seq_along(c)) {
                    subr[[i]] <- list(sample(nrow(test.mat), n, replace=TRUE),
                                      sample(all.values, n, replace=TRUE))
                    to.use <- !duplicated(subr[[i]][[1]], fromLast=TRUE) # Last elements overwrite earlier elements.
                    ref[subr[[i]][[1]][to.use],c[i]] <- subr[[i]][[2]][to.use]
                }
            }

            out <- .Call(cxxfun, test.mat, it, c, subr)
            if (is.matrix(test.mat)) { 
                testthat::expect_true(is.matrix(out))
            } else {
                testthat::expect_s4_class(out, class(test.mat))
                out <- as.matrix(out)
                dimnames(out) <- NULL
            }

            testthat::expect_identical(ref, out)
        }
    }
    return(invisible(NULL))
}

check_integer_output_indexed <- function(FUN, ..., N) {
    .check_output_indexed(FUN=FUN, ..., N=N,cxxfun=cxx_test_integer_output_indexed, fill=0L)
} 

check_logical_output_indexed <- function(FUN, ..., N) {
    .check_output_indexed(FUN=FUN, ..., N=N,cxxfun=cxx_test_logical_output_indexed, fill=FALSE)
} 

check_numeric_output_indexed <- function(FUN, ..., N) {
    .check_output_indexed(FUN=FUN, ..., N=N,cxxfun=cxx_test_numeric_output_indexed, fill=0)
} 

check_character_output_indexed <- function(FUN, ..., N) {
    .check_output_indexed(FUN=FUN, ..., N=N,cxxfun=cxx_test_character_output_indexed, fill="")
} 

###############################

.check_execution_order <- function(FUN, cxxfun, type) {
    # Checking that the output function in '.Call' does not overwrite the 
    # underlying HDF5 file and change the values of other HDF5Matrix objects. 
    mats <- list(FUN(), FUN(), FUN())
    ref.mats <- lapply(mats, as.matrix)
    new.mats <- lapply(mats, function(x) .Call(cxxfun, x, 1L, NULL)[[1]]) 

    for (i in seq_along(mats)) { 
        out <- new.mats[[i]]
        testthat::expect_s4_class(out, "HDF5Matrix")
        ref <- ref.mats[[i]]
        testthat::expect_identical(ref, as.matrix(out))
        original <- mats[[i]]
        testthat::expect_identical(ref, as.matrix(original))
        testthat::expect_true(original@seed@filepath!=out@seed@filepath)
    }

    # Checking that the old and realized files are in the log.
    testthat::expect_message(log <- HDF5Array::showHDF5DumpLog())
    for (mode in c(TRUE, FALSE)) {
        for (i in seq_along(mats)) {
            if (mode) {
                current <- new.mats[[i]]
            } else {
                current <- mats[[i]]
            }   

            j <- which(log$name==current@seed@name & log$file==current@seed@filepath)
            testthat::expect_true(length(j)==1L)
            testthat::expect_identical(type, log$type[j])
            testthat::expect_identical(sprintf("%ix%i", nrow(current), ncol(current)), log$dims[j])
        }
    }
    return(invisible(NULL))
}

check_integer_order <- function(FUN, cxxfun) {
    .check_execution_order(FUN, cxxfun=cxx_test_integer_output, type="integer")
}

check_numeric_order <- function(FUN, cxxfun) {
    .check_execution_order(FUN, cxxfun=cxx_test_numeric_output, type="double")
}

check_logical_order <- function(FUN, cxxfun) {
    .check_execution_order(FUN, cxxfun=cxx_test_logical_output, type="logical")
}

check_character_order <- function(FUN, cxxfun) {
    .check_execution_order(FUN, cxxfun=cxx_test_character_output, type="character")
}

###############################

.check_converted_output <- function(FUN, ..., cxxfun, rfun) { 
    for (i in 1:3) { 
        test.mat <- FUN(...)
        
        out <- .Call(cxxfun, test.mat, i)
        if (is.matrix(test.mat)) { 
            testthat::expect_true(is.matrix(out[[1]]))
        } else {
            testthat::expect_s4_class(out[[1]], class(out[[1]]))
            out[[1]] <- as.matrix(out[[1]])
        }
        ref <- as.matrix(test.mat)
        testthat::expect_identical(rfun(ref), out[[1]])
            
        # Checking that the converted getters work properly too.
        ref2 <- matrix(out[[2]], nrow(ref), ncol(ref), byrow=(i==2L))
        testthat::expect_identical(ref, ref2)
    }
    return(invisible(NULL))
}

check_numeric_converted_output  <- function(FUN, ...) {
    .check_converted_output(FUN=FUN, ..., cxxfun=cxx_test_numeric_to_integer_output, 
                            rfun=function(x) {
                                storage.mode(x) <- "integer" 
                                return(x)
                            })
}

check_integer_converted_output <- function(FUN, ...) {
    .check_converted_output(FUN=FUN, ..., cxxfun=cxx_test_integer_to_numeric_output, 
                            rfun=function(x) { 
                                storage.mode(x) <- "double"
                                return(x)
                            })
}

check_logical_converted_output <- function(FUN, ...) {
    .check_converted_output(FUN=FUN, ..., cxxfun=cxx_test_logical_to_numeric_output, 
                            rfun=function(x) { 
                                storage.mode(x) <- "double"
                                return(x)
                            })
    .check_converted_output(FUN=FUN, ..., cxxfun=cxx_test_logical_to_integer_output, 
                            rfun=function(x) { 
                                storage.mode(x) <- "integer"
                                return(x)
                            })
}

###############################

check_output_mode <- function(incoming, ..., simplify, preserve.zero) {
    all.modes <- .Call(cxx_get_all_modes)
    if (is.character(incoming)) { 
        out <- .Call(cxx_select_output_by_mode, incoming, simplify, preserve.zero)
    } else {
        out <- .Call(cxx_select_output_by_sexp, incoming(...), simplify, preserve.zero)
    }
    return(names(all.modes)[all.modes==out])
}

###############################

.check_edge_output_errors <- function(x, cxxfun) {
    testthat::expect_true(.Call(cxxfun, x, 0L))

    testthat::expect_error(.Call(cxxfun, x, 1L), "row index out of range")
    testthat::expect_error(.Call(cxxfun, x, -1L), "column index out of range")        
    testthat::expect_error(.Call(cxxfun, x, 2L), "column start index is greater than column end index")
    testthat::expect_error(.Call(cxxfun, x, -2L), "row start index is greater than row end index")
    testthat::expect_error(.Call(cxxfun, x, 3L), "column end index out of range")
    testthat::expect_error(.Call(cxxfun, x, -3L), "row end index out of range")
    
    testthat::expect_error(.Call(cxxfun, x, 4L), "row index out of range")
    testthat::expect_error(.Call(cxxfun, x, -4L), "column index out of range")        
    testthat::expect_error(.Call(cxxfun, x, 5L), "column start index is greater than column end index")
    testthat::expect_error(.Call(cxxfun, x, -5L), "row start index is greater than row end index")
    testthat::expect_error(.Call(cxxfun, x, 6L), "column end index out of range")
    testthat::expect_error(.Call(cxxfun, x, -6L), "row end index out of range")

    return(invisible(NULL))
}

check_integer_edge_output_errors <- function(FUN, ...) {
    .check_edge_output_errors(FUN(...), cxxfun=cxx_test_integer_edge_output)
}

check_logical_edge_output_errors <- function(FUN, ...) {
    .check_edge_output_errors(FUN(...), cxxfun=cxx_test_logical_edge_output)
}

check_numeric_edge_output_errors <- function(FUN, ...) {
    .check_edge_output_errors(FUN(...), cxxfun=cxx_test_numeric_edge_output)
}

check_character_edge_output_errors <- function(FUN, ...) {
    .check_edge_output_errors(FUN(...), cxxfun=cxx_test_character_edge_output)
}

