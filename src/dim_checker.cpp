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
    if (d[0]<0 || d[1]<0) { 
        throw std::runtime_error("dimensions should be non-negative"); 
    }
    nrow=d[0];
    ncol=d[1];
    return;
}

void dim_checker::check_dimension(size_t i, size_t dim, const char* msg) {
    if (i >= dim) {
        std::stringstream err;
        err << msg << " index out of range";
        throw std::runtime_error(err.str().c_str());
    }
    return;
}

void dim_checker::check_rowargs(size_t r) const {
    dim_checker::check_dimension(r, nrow, "row");
    return;
}

void dim_checker::check_colargs(size_t c) const {
    dim_checker::check_dimension(c, ncol, "column");
    return;
}

void dim_checker::check_subset(size_t first, size_t last, size_t dim, const char* msg) {
    if (last < first) {
        std::stringstream err;
        err << msg << " start index is greater than " << msg << " end index";
        throw std::runtime_error(err.str().c_str());
    } else if (last > dim) {
        std::stringstream err;
        err << msg << " end index out of range";
        throw std::runtime_error(err.str().c_str());
    }
    return;    
}

void dim_checker::check_rowargs(size_t r, size_t first, size_t last) const {
    check_rowargs(r);
    dim_checker::check_subset(first, last, ncol, "column");
    return;
}


void dim_checker::check_colargs(size_t c, size_t first, size_t last) const {
    check_colargs(c);
    dim_checker::check_subset(first, last, nrow, "row");
    return;
}

void dim_checker::check_oneargs(size_t r, size_t c) const {
    check_rowargs(r);
    check_colargs(c);
    return;
}

void check_indices(Rcpp::IntegerVector::iterator it, size_t n, size_t dim, const char * msg) { 
    if (n==0) { return; }

    int last=*(it++);
    for (size_t i=1; i<n; ++i, ++it) {
        dim_checker::check_dimension(*it, dim, msg);
        if (*it <= last) {
            std::stringstream err;
            err << msg << " indices are not strictly increasing";
            throw std::runtime_error(err.str().c_str());
        }
    }

    return;
}

void dim_checker::check_row_indices(Rcpp::IntegerVector::iterator it, size_t n) {
    check_indices(it, n, nrow, "row");
    return;
}

void dim_checker::check_col_indices(Rcpp::IntegerVector::iterator it, size_t n) {
    check_indices(it, n, ncol, "column");
    return;
}

}
