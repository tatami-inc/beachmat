#ifndef BEACHMAT_ORDINARY_READER_H
#define BEACHMAT_ORDINARY_READER_H

#include "Rcpp.h"
#include "dim_checker.h"
#include "utils.h"

namespace beachmat {

template <class V>
class ordinary_reader : public dim_checker {
public:
    ~ordinary_reader() = default;
    ordinary_reader(const ordinary_reader&) = default;
    ordinary_reader& operator=(const ordinary_reader&) = default;
    ordinary_reader(ordinary_reader&&) = default;
    ordinary_reader& operator=(ordinary_reader&&) = default;

    ordinary_reader(Rcpp::RObject input) {
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

    typename V::const_iterator get_col(size_t c, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        return mat.begin() + (c * this->nrow) + first;
    }

    template <class Iter>
    void get_row(size_t r, Iter work, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        auto src=mat.begin() + first * (this->nrow) + r;
        for (size_t col=first; col<last; ++col, src+=(this->nrow), ++work) { (*work)=(*src); }
        return;
    }
private:
    V mat;
};

}

#endif
