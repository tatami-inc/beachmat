#ifndef BEACHMAT_LOGICAL_MATRIX_H
#define BEACHMAT_LOGICAL_MATRIX_H

#include "input/LIN_matrix.h"
#include "output/LIN_output.h"

#include <memory>
#include <stdexcept>

namespace beachmat {

/*********
 * INPUT *
 *********/

/* Virtual base class for logical matrices. */

typedef lin_matrix<int, Rcpp::LogicalVector> logical_matrix;

std::unique_ptr<logical_matrix> create_logical_matrix_internal(const Rcpp::RObject&, bool); 

/* Simple logical matrix */

typedef simple_lin_matrix<int, Rcpp::LogicalVector> simple_logical_matrix;

/* lgeMatrix */

typedef dense_lin_matrix<int, Rcpp::LogicalVector> dense_logical_matrix;

/* lgCMatrix */

template<>
int Csparse_reader<int, Rcpp::LogicalVector>::get_empty() { return 0; }

typedef Csparse_lin_matrix<int, Rcpp::LogicalVector> Csparse_logical_matrix;

/* DelayedMatrix */

template<>
std::unique_ptr<logical_matrix> delayed_lin_reader<int, Rcpp::LogicalVector>::generate_seed(Rcpp::RObject incoming) {
    return create_logical_matrix_internal(incoming, false);
}

typedef delayed_lin_matrix<int, Rcpp::LogicalVector> delayed_logical_matrix;

/* Unknown matrix */

typedef unknown_lin_matrix<int, Rcpp::LogicalVector> unknown_logical_matrix;

/* External matrix */

typedef external_lin_matrix<int, Rcpp::LogicalVector> external_logical_matrix;

/* Dispatcher */

inline std::unique_ptr<logical_matrix> create_logical_matrix_internal(const Rcpp::RObject& incoming, bool delayed) { 
    if (incoming.isS4()) {
        std::string ctype=get_class(incoming);
        if (ctype=="lgeMatrix") { 
            return std::unique_ptr<logical_matrix>(new dense_logical_matrix(incoming));
        } else if (ctype=="lgCMatrix") { 
            return std::unique_ptr<logical_matrix>(new Csparse_logical_matrix(incoming));
        } else if (delayed && ctype=="DelayedMatrix") { 
            return std::unique_ptr<logical_matrix>(new delayed_logical_matrix(incoming));
        } else if (has_external_support(incoming)) {
            return std::unique_ptr<logical_matrix>(new external_logical_matrix(incoming));
        }
        return std::unique_ptr<logical_matrix>(new unknown_logical_matrix(incoming));
    } 
    return std::unique_ptr<logical_matrix>(new simple_logical_matrix(incoming));
}

inline std::unique_ptr<logical_matrix> create_logical_matrix(const Rcpp::RObject& incoming) {
    return create_logical_matrix_internal(incoming, true);
}

/**********
 * OUTPUT *
 **********/

/* Virtual base class for output logical matrices. */

typedef lin_output<int, Rcpp::LogicalVector> logical_output;

/* Simple output logical matrix */

typedef simple_lin_output<int, Rcpp::LogicalVector> simple_logical_output;

/* Sparse output logical matrix */

template<>
int Csparse_writer<int, Rcpp::LogicalVector>::get_empty() { return 0; }

typedef sparse_lin_output<int, Rcpp::LogicalVector> sparse_logical_output;

/* Output dispatchers */

inline std::unique_ptr<logical_output> create_logical_output(int nrow, int ncol, const output_param& param) {
    switch (param.get_mode()) {
        case SIMPLE:
            return std::unique_ptr<logical_output>(new simple_logical_output(nrow, ncol));
        case SPARSE:
            return std::unique_ptr<logical_output>(new sparse_logical_output(nrow, ncol));
        default:
            throw std::runtime_error("unsupported output mode for logical matrices");
    }
}

}

#endif
