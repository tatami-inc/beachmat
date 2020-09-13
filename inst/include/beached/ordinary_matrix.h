#ifndef BEACHMAT_ORDINARY_MATRIX_H
#define BEACHMAT_ORDINARY_MATRIX_H

#include "Rcpp.h"
#include "any_matrix.h"
#include "utils.h"

namespace beachmat {

template <class V, typename T = typename V::stored_type>
class ordinary_matrix : public any_matrix<T> {
public:
    virtual ~ordinary_matrix() = default;
    ordinary_matrix(const ordinary_matrix&) = default;
    ordinary_matrix& operator=(const ordinary_matrix&) = default;
    ordinary_matrix(ordinary_matrix&&) = default;
    ordinary_matrix& operator=(ordinary_matrix&&) = default;

    ordinary_matrix(Rcpp::RObject input) {
        if (!input.hasAttribute("dim")) { 
            throw std::runtime_error("matrix object should have 'dim' attribute"); 
        }
        this->fill_dims(input.attr("dim"));
        const size_t& NC=this->ncol; 

        if (input.sexp_type()!=mat.sexp_type()) { 
            throw std::runtime_error(std::string("matrix should be ") + translate_type(mat.sexp_type()));
        }
        mat=input;
        if (static_cast<size_t>(mat.size())!=(this->nrow)*NC) {
            throw std::runtime_error("length of matrix is inconsistent with its dimensions"); 
        }
        return;
    }

    T* get_col(size_t c, T* work, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        return mat.begin() + (c * this->nrow) + first;
    }

    T* get_row(size_t r, T* work, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        auto src=mat.begin() + first * (this->nrow) + r;
        auto copy=work;
        for (size_t col=first; col<last; ++col, src+=(this->nrow), ++work) { (*work)=(*src); }
        return copy;
    }
private:
    V mat;
};

}

#endif
