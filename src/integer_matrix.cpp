#include "integer_matrix.h"

namespace beachmat {

/* DelayedMatrix input methods. */

template<>
std::unique_ptr<integer_matrix> delayed_lin_matrix<int, Rcpp::IntegerVector>::generate_seed(Rcpp::RObject incoming) {
    // incoming should be a "DelayedMatrix" object, not the seed within it!
    bool isokay=false;
    Rcpp::RObject seed(get_safe_slot(incoming, "seed"));

    if (seed.isS4()) { 
        std::string ctype=get_class(seed);
        if (ctype=="RleMatrix") {
            isokay=true;
            incoming=seed;
        } else if (ctype=="HDF5ArraySeed") {
            isokay=true;
            incoming=delayed_seed_to_HDF5Matrix(seed);
        }
    } else {
        isokay=true;
        incoming=seed;
    }

    if (isokay) {
        return create_integer_matrix(incoming);
    } else {
        return nullptr;
    }
} 
    
/* HDF5 integer output methods. */

template<>
int HDF5_output<int, INTSXP>::get_empty() const { return 0; }

template<>
Rcpp::RObject HDF5_output<int, INTSXP>::get_firstval() { 
    int first;
    extract_one(0, 0, &first);
    return Rcpp::IntegerVector::create(first);
}

/* Dispatch definition */

std::unique_ptr<integer_matrix> create_integer_matrix(const Rcpp::RObject& incoming) { 
    if (incoming.isS4()) { 
        std::string ctype=get_class(incoming);
        if (ctype=="HDF5Matrix") { 
            return std::unique_ptr<integer_matrix>(new HDF5_integer_matrix(incoming));
        } else if (ctype=="RleMatrix") {
            return std::unique_ptr<integer_matrix>(new Rle_integer_matrix(incoming));
        } else if (ctype=="DelayedMatrix") {
            return std::unique_ptr<integer_matrix>(new delayed_integer_matrix(incoming));  
        }
        std::stringstream err;
        err << "unsupported class '" << ctype << "' for integer_matrix";
        throw std::runtime_error(err.str().c_str());
    } 
    return std::unique_ptr<integer_matrix>(new simple_integer_matrix(incoming));
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
