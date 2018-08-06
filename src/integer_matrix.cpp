#include "integer_matrix.h"

namespace beachmat {

/* HDF5Matrix input methods. */

template<>
int HDF5_lin_matrix<int, Rcpp::IntegerVector, INTSXP>::get(size_t r, size_t c) {
    int out;
    reader.extract_one(r, c, &out, H5::PredType::NATIVE_INT32);
    return out; 
}

/* DelayedMatrix input methods. */

std::unique_ptr<integer_matrix> create_integer_matrix_internal(const Rcpp::RObject&, bool);

template<>
std::unique_ptr<integer_matrix> delayed_lin_reader<int, Rcpp::IntegerVector>::generate_seed(Rcpp::RObject incoming) {
    return create_integer_matrix_internal(incoming, false);
} 

/* HDF5 integer output methods. */

template<>
int HDF5_writer<int, INTSXP>::get_empty() { return 0; }

template<>
Rcpp::RObject HDF5_writer<int, INTSXP>::get_firstval() { 
    int first;
    extract_one(0, 0, &first);
    return Rcpp::IntegerVector::create(first);
}

/* Dispatch definition */

std::unique_ptr<integer_matrix> create_integer_matrix_internal(const Rcpp::RObject& incoming, bool delayed) { 
    if (incoming.isS4()) { 
        std::string ctype=get_class(incoming);
        if (ctype=="HDF5Matrix") { 
            return std::unique_ptr<integer_matrix>(new HDF5_integer_matrix(incoming));
        } else if (delayed && ctype=="DelayedMatrix") {
            return std::unique_ptr<integer_matrix>(new delayed_integer_matrix(incoming));  
        }
        return std::unique_ptr<integer_matrix>(new unknown_integer_matrix(incoming));
    } 
    return std::unique_ptr<integer_matrix>(new simple_integer_matrix(incoming));
}

std::unique_ptr<integer_matrix> create_integer_matrix(const Rcpp::RObject& incoming) { 
    return create_integer_matrix_internal(incoming, true);
}

/* Output dispatch definition */

std::unique_ptr<integer_output> create_integer_output(int nrow, int ncol, const output_param& param) {
    switch (param.get_mode()) {
        case SIMPLE:
            return std::unique_ptr<integer_output>(new simple_integer_output(nrow, ncol));
        case HDF5:
            return std::unique_ptr<integer_output>(new HDF5_integer_output(nrow, ncol,
                        param.get_chunk_nrow(), param.get_chunk_ncol(), param.get_compression()));
        default:
            throw std::runtime_error("unsupported output mode for integer matrices");
    }
}

}
