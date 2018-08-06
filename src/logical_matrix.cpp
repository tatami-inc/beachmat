#include "logical_matrix.h"

namespace beachmat {
    
/* Csparse logical input methods. */

template<>
int Csparse_reader<int, Rcpp::LogicalVector>::get_empty() { return 0; }

/* HDF5Matrix input methods. */

template<>
int HDF5_lin_reader<int, Rcpp::LogicalVector, LGLSXP>::get(size_t r, size_t c) {
    int out;
    mat.extract_one(r, c, &out, H5::PredType::NATIVE_INT32);
    return out; 
}

/* DelayedMatrix input methods. */

std::unique_ptr<logical_matrix> create_logical_matrix_internal(const Rcpp::RObject&, bool); 

template<>
std::unique_ptr<logical_matrix> delayed_lin_reader<int, Rcpp::LogicalVector>::generate_seed(Rcpp::RObject incoming) {
    return create_logical_matrix_internal(incoming, false);
} 

/* Sparse logical output methods. */

template<>
int Csparse_writer<int, Rcpp::LogicalVector>::get_empty() { return 0; }

/* HDF5 logical output methods. */

template<>
int HDF5_output<int, LGLSXP>::get_empty() { return 0; }

template<>
Rcpp::RObject HDF5_output<int, LGLSXP>::get_firstval() { 
    int first;
    extract_one(0, 0, &first);
    return Rcpp::LogicalVector::create(first);
}

/* Dispatch definition */

std::unique_ptr<logical_matrix> create_logical_matrix_internal(const Rcpp::RObject& incoming, bool delayed) { 
    if (incoming.isS4()) {
        std::string ctype=get_class(incoming);
        if (ctype=="lgeMatrix") { 
            return std::unique_ptr<logical_matrix>(new dense_logical_matrix(incoming));
        } else if (ctype=="lgCMatrix") { 
            return std::unique_ptr<logical_matrix>(new Csparse_logical_matrix(incoming));
        } else if (ctype=="lgTMatrix") {
            throw std::runtime_error("lgTMatrix not supported, convert to lgCMatrix");
        } else if (ctype=="HDF5Matrix") {
            return std::unique_ptr<logical_matrix>(new HDF5_logical_matrix(incoming));
        } else if (delayed && ctype=="DelayedMatrix") { 
            return std::unique_ptr<logical_matrix>(new delayed_logical_matrix(incoming));
        }
        return std::unique_ptr<logical_matrix>(new unknown_logical_matrix(incoming));
    } 
    return std::unique_ptr<logical_matrix>(new simple_logical_matrix(incoming));
}

std::unique_ptr<logical_matrix> create_logical_matrix(const Rcpp::RObject& incoming) {
    return create_logical_matrix_internal(incoming, true);
}

/* Output dispatch definition */

std::unique_ptr<logical_output> create_logical_output(int nrow, int ncol, const output_param& param) {
    switch (param.get_mode()) {
        case SIMPLE:
            return std::unique_ptr<logical_output>(new simple_logical_output(nrow, ncol));
        case SPARSE:
            return std::unique_ptr<logical_output>(new sparse_logical_output(nrow, ncol));
        case HDF5:
            return std::unique_ptr<logical_output>(new HDF5_logical_output(nrow, ncol,
                        param.get_chunk_nrow(), param.get_chunk_ncol(), param.get_compression()));
        default:
            throw std::runtime_error("unsupported output mode for logical matrices");
    }
}

}
