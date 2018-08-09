#include "beachtest.h"

extern "C" {

SEXP set_class_by_sexp(SEXP incoming, SEXP simplify, SEXP preserve_zero) {
    beachmat::output_param op(incoming, Rf_asLogical(simplify), Rf_asLogical(preserve_zero));

    Rcpp::StringVector out(1);
    switch (op.get_mode()) {
        case beachmat::SIMPLE:
            out[0]="simple";
            return out;
        case beachmat::SPARSE:
            out[0]="sparse";
            return out;
        case beachmat::HDF5:
            out[0]="HDF5";
            return out;
    }
    throw std::runtime_error("invalid output type");
}

}

