#ifndef BEACHMAT_OUTPUT_PARAM_H
#define BEACHMAT_OUTPUT_PARAM_H

#include "Rcpp.h"

#include "../utils/utils.h"

namespace beachmat {

class output_param {
public:
    output_param (matrix_type m) : mode(m) {}

    output_param (const Rcpp::RObject& in, bool simplify, bool preserve_zero) : 
        output_param(robject_to_matrix_class(in), simplify, preserve_zero) {}

    output_param(matrix_type m, bool simplify, bool preserve_zero=false) : output_param(m) {
        switch (mode) {
            case SIMPLE:
                break;
            case DENSE:
                mode=SIMPLE;
                break;
            case SPARSE:
                if (preserve_zero) { break; } // keeping sparse, if preserve_zero is true.
            default:
                mode=SIMPLE;
        }
        return; 
    }

    matrix_type get_mode() const {
        return mode;
    }    
private:
    matrix_type mode;

    static matrix_type robject_to_matrix_class(const Rcpp::RObject& in) {
        if (in.isS4()) {
            auto curclass=get_class(in);
            if (curclass=="DelayedMatrix") { 
                return DELAYED;
            } else if (!curclass.empty() && curclass.substr(1)=="gCMatrix") {
                return SPARSE;
            } else if (!curclass.empty() && curclass.substr(1)=="geMatrix") {
                return DENSE;
            }
            return UNKNOWN;
        }
        return SIMPLE;
    }
};

}

#endif
