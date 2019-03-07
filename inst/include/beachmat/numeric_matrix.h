#ifndef BEACHMAT_NUMERIC_MATRIX_H
#define BEACHMAT_NUMERIC_MATRIX_H

#include "LIN_matrix.h"
#include "LIN_output.h"

#include <sstream>
#include <stdexcept>

namespace beachmat { 

/*********
 * INPUT *
 *********/ 

std::unique_ptr<numeric_matrix> create_numeric_matrix_internal(const Rcpp::RObject&, bool); 

/* Virtual base class for numeric matrices. */

typedef lin_matrix<double, Rcpp::NumericVector> numeric_matrix;

/* Simple numeric matrix */

typedef simple_lin_matrix<double, Rcpp::NumericVector> simple_numeric_matrix;

/* dgeMatrix */

typedef dense_lin_matrix<double, Rcpp::NumericVector> dense_numeric_matrix;

/* dgCMatrix */

template<>
double Csparse_reader<double, Rcpp::NumericVector>::get_empty() { return 0; }

typedef Csparse_lin_matrix<double, Rcpp::NumericVector> Csparse_numeric_matrix;

/* DelayedMatrix */

template<>
std::unique_ptr<numeric_matrix> delayed_lin_reader<double, Rcpp::NumericVector>::generate_seed(Rcpp::RObject incoming) {
    return create_numeric_matrix_internal(incoming, false);
}

typedef delayed_lin_matrix<double, Rcpp::NumericVector> delayed_numeric_matrix;

/* Unknown matrix */

typedef unknown_lin_matrix<double, Rcpp::NumericVector> unknown_numeric_matrix;

/* External matrix */

typedef external_lin_matrix<double, Rcpp::NumericVector> external_numeric_matrix;

/* Dispatcher */

inline std::unique_ptr<numeric_matrix> create_numeric_matrix_internal(const Rcpp::RObject& incoming, bool delayed) { 
    if (incoming.isS4()) {
        std::string ctype=get_class(incoming);
        if (ctype=="dgeMatrix") { 
            return std::unique_ptr<numeric_matrix>(new dense_numeric_matrix(incoming));
        } else if (ctype=="dgCMatrix") { 
            return std::unique_ptr<numeric_matrix>(new Csparse_numeric_matrix(incoming));
        } else if (delayed && ctype=="DelayedMatrix") { 
            return std::unique_ptr<numeric_matrix>(new delayed_numeric_matrix(incoming));
        } else if (has_external_support(incoming)) {
            return std::unique_ptr<numeric_matrix>(new external_numeric_matrix(incoming));
        }
        return std::unique_ptr<numeric_matrix>(new unknown_numeric_matrix(incoming));
    } 
    return std::unique_ptr<numeric_matrix>(new simple_numeric_matrix(incoming));
}

inline std::unique_ptr<numeric_matrix> create_numeric_matrix(const Rcpp::RObject& incoming) { 
    return create_numeric_matrix_internal(incoming, true);
}

/**********
 * OUTPUT *
 **********/ 

/* Virtual base class for output numeric matrices. */

typedef lin_output<double, Rcpp::NumericVector> numeric_output;

/* Simple output numeric matrix */

typedef simple_lin_output<double, Rcpp::NumericVector> simple_numeric_output;

/* Sparse output numeric matrix */

template<>
double Csparse_writer<double, Rcpp::NumericVector>::get_empty() { return 0; }

typedef sparse_lin_output<double, Rcpp::NumericVector> sparse_numeric_output;

/* Output dispatchers */

inline std::unique_ptr<numeric_output> create_numeric_output(int nrow, int ncol, const output_param& param) {
    switch (param.get_mode()) {
        case SIMPLE:
            return std::unique_ptr<numeric_output>(new simple_numeric_output(nrow, ncol));
        case SPARSE:
            return std::unique_ptr<numeric_output>(new sparse_numeric_output(nrow, ncol));
        default:
            throw std::runtime_error("unsupported output mode for numeric matrices");
    }
}

}

#endif
