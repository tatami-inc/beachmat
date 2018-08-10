#include "beachtest.h"

extern "C" {

SEXP set_class_by_sexp(SEXP incoming, SEXP simplify, SEXP preserve_zero) {
    BEGIN_RCPP
    beachmat::output_param op(incoming, Rf_asLogical(simplify), Rf_asLogical(preserve_zero));
    return translate_class(op.get_mode());
    END_RCPP
}

}

