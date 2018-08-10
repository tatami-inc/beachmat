#include "beachtest.h"

Rcpp::StringVector translate_class(beachmat::matrix_type X) {
    Rcpp::StringVector out(1);
    switch (X) {
        case beachmat::SIMPLE:
            out[0]="simple";
            return out;
        case beachmat::DENSE:
            out[0]="dense";
            return out;
        case beachmat::SPARSE:
            out[0]="sparse";
            return out;
        case beachmat::UNKNOWN:
            out[0]="unknown";
            return out;
        case beachmat::DELAYED:
            out[0]="delayed";
            return out;
        case beachmat::HDF5:
            out[0]="HDF5";
            return out;
    }
    throw std::runtime_error("invalid output type");
}

extern "C" {

SEXP get_class_integer(SEXP incoming) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(incoming);
    return translate_class(ptr->get_matrix_type());	
    END_RCPP
}

SEXP get_class_numeric(SEXP incoming) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(incoming);
    return translate_class(ptr->get_matrix_type());	
    END_RCPP
}

SEXP get_class_logical(SEXP incoming) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(incoming);
    return translate_class(ptr->get_matrix_type());	
    END_RCPP
}

SEXP get_class_character(SEXP incoming) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(incoming);
    return translate_class(ptr->get_matrix_type());	
    END_RCPP
}

}

