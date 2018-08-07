spawn_row_ordering <- function(NROW) {
    rranges <- list(forward=seq_len(NROW), random=sample(NROW), subset=sample(NROW, NROW/2L))
    rranges$reverse <- rev(rranges$forward)
    rranges
}

spawn_col_ordering <- function(NCOL) {
    cranges <- list(forward=seq_len(NCOL), random=sample(NCOL), subset=sample(NCOL, NCOL/2L))
    cranges$reverse <- rev(cranges$forward)
    cranges
}
