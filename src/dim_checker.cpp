#include "dim_checker.h"

namespace beachmat { 

dim_checker::dim_checker() : nrow(0), ncol(0) {}

dim_checker::dim_checker(size_t nr, size_t nc) : nrow(nr), ncol(nc) {}

dim_checker::~dim_checker() {}

size_t dim_checker::get_nrow() const { return nrow; }

size_t dim_checker::get_ncol() const { return ncol; }

void dim_checker::fill_dims(const Rcpp::RObject& dims) {
    Rcpp::IntegerVector d;
    if (dims.sexp_type()!=d.sexp_type() || (d=dims).size()!=2) { 
        throw std::runtime_error("matrix dimensions should be an integer vector of length 2");
    }
    if (d[0]<0 || d[1]<0) { throw std::runtime_error("dimensions should be non-negative"); }
    nrow=d[0];
    ncol=d[1];
    return;
}

void dim_checker::check_rowargs(size_t r) const {
    if (r>=nrow) {
        throw std::runtime_error("row index out of range");
    }
    return;
}

void dim_checker::check_rowargs(size_t r, size_t first, size_t last) const {
    check_rowargs(r);
    if (last < first) {
        throw std::runtime_error("column start index is greater than column end index");
    } else if (last > ncol) {
        throw std::runtime_error("column end index out of range");
    }
    return;    
}

void dim_checker::check_colargs(size_t c) const {
    if (c>=ncol) {
        throw std::runtime_error("column index out of range");
    }
    return;
}

void dim_checker::check_colargs(size_t c, size_t first, size_t last) const {
    check_colargs(c);
    if (last < first) {
        throw std::runtime_error("row start index is greater than row end index");
    } else if (last > nrow) {
        throw std::runtime_error("row end index out of range");
    }
    return;
}

void dim_checker::check_oneargs(size_t r, size_t c) const {
    if (c>=ncol || r>=nrow) {
        throw std::runtime_error("column or row indices out of range");
    }
    return;
}

}
