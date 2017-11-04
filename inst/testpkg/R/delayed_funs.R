delayed_funs <- function(basefun) 
# Creates a list of supported DelayedMatrix operations 
# that are handled natively in beachmat, without the 
# need for explicit realization.
{
    # Subsetting:
    .choose_half <- function(n) {
        sample(n, ceiling(n/2))
    }

    subset_row_fun <- function(...) {
        out <- DelayedArray::DelayedArray(basefun(...))
        out[.choose_half(nrow(out)),]
    }

    subset_col_fun <- function(...) {
        out <- DelayedArray::DelayedArray(basefun(...))
        out[,.choose_half(ncol(out))]
    }

    subset_both_fun <- function(...) {
        out <- DelayedArray::DelayedArray(basefun(...))
        out[.choose_half(nrow(out)),.choose_half(ncol(out))]
    }

    # Dimnaming
    .namers <- function(n, pre) {
        paste0(pre, seq_len(n))
    }

    name_row_fun <- function(...) { 
        out <- DelayedArray::DelayedArray(basefun(...))
        rownames(out) <- .namers(nrow(out), "Gene")
        return(out)
    }

    name_col_fun <- function(...) { 
        out <- DelayedArray::DelayedArray(basefun(...))
        colnames(out) <- .namers(ncol(out), "Cell")
        return(out)
    }

    name_both_fun <- function(...) {
        out <- DelayedArray::DelayedArray(basefun(...))
        rownames(out) <- .namers(nrow(out), "Gene")
        colnames(out) <- .namers(ncol(out), "Cell")
        return(out)
    }

    # Transposition
    trans_fun <- function(...) { 
        out <- DelayedArray::DelayedArray(basefun(...))
        DelayedArray::t(out)
    }

    trans_subset_fun <- function(...) {
        DelayedArray::t(subset_both_fun(...))
    }

    return(list(SR=subset_row_fun,
                SC=subset_col_fun,
                SB=subset_both_fun,
                NR=name_row_fun,
                NC=name_col_fun,
                NB=name_both_fun,
                TR=trans_fun,
                TRS=trans_subset_fun))
}
