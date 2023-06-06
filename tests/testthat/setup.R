am_i_ok <- function(ref, ptr, exact=TRUE) {
    expect_identical(dim(ref), beachmat:::tatami_dim(ptr))
    test <- if (exact) expect_identical else expect_equal
    for (i in seq_len(ncol(ref))) {
        expected <- ref[,i]
        if (!is.double(expected)) { 
            expected <- as.double(expected) 
        }
        test(expected, beachmat:::tatami_column(ptr, i))
    }

    # Checking for thread-correct processing.
    expect_equal(beachmat:::tatami_row_sums(ptr, 2), Matrix::rowSums(ref))
    expect_equal(beachmat:::tatami_column_sums(ptr, 2), Matrix::colSums(ref))
}
