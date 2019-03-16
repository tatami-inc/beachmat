#ifndef BEACHTEST_GET_CONST_H
#define BEACHTEST_GET_CONST_H
#include "beachtest.h"

template <class T, class O, class M>  // M is automatically deduced.
O get_const_all (M ptr, Rcpp::IntegerVector ordering) {
    const size_t nrows=ptr->get_nrow();
    O output(nrows, ordering.size());

    const bool use_const=(ptr->col_raw_type()=="dense");
    auto raws=ptr->set_up_raw();
    T alternative(nrows);
    auto it=alternative.begin();
    size_t c=0;

    for (auto o : ordering) {
        if (use_const) {
            ptr->get_col_raw(o-1, raws);
            it=raws.get_values_start();
        } else {
            ptr->get_col(o-1, it);
        }

        auto outcol=output.column(c);
        std::copy(it, it + nrows, outcol.begin());
        ++c;
    }
    return output;
}

template <class T, class O, class M>  
O get_const_slice (M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerVector rows) {
    if (rows.size()!=2) { 
        throw std::runtime_error("'rows' should be an integer vector of length 2"); 
    }
    const int rstart=rows[0]-1, rend=rows[1];
    const int nrows=rend-rstart;    
    O output(nrows, ordering.size());

    const bool use_const=(ptr->col_raw_type()=="dense");
    auto raws=ptr->set_up_raw();
    T alternative(nrows);
    auto it=alternative.begin();
    size_t c=0;

    for (auto o : ordering) {
        if (use_const) {
            ptr->get_col_raw(o-1, raws, rstart, rend);
            it=raws.get_values_start();
        } else {
            ptr->get_col(o-1, it, rstart, rend);
        }

        auto curcol=output.column(c);
        std::copy(it, it+nrows, curcol.begin());
        ++c;
    }

    return output;
}

template <class T, class M>  
Rcpp::List get_const_varslice (M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerMatrix rows) {
    if (rows.ncol()!=2) { 
        throw std::runtime_error("'rows' should be an integer matrix with two columns"); 
    }
    if (rows.nrow()!=ordering.size()) {
        throw std::runtime_error("'nrow(rows)' should be equal to 'length(ordering)'");
    }
    Rcpp::List output(ordering.size());

    const bool use_const=(ptr->col_raw_type()=="dense");
    auto raws=ptr->set_up_raw();
    T alternative(ptr->get_nrow());
    auto it=alternative.begin();
    size_t c=0;

    for (auto o : ordering) {
        auto cur_bounds=rows.row(c);
        int left=cur_bounds[0]-1, right=cur_bounds[1];

        if (use_const) {
            ptr->get_col_raw(o-1, raws, left, right);
            it=raws.get_values_start();
        } else {
            ptr->get_col(o-1, it, left, right);
        }

        T out(right-left);
        std::copy(it, it+right-left, out.begin());
        output[c]=out;
        ++c;
    }

    return output;
}

#endif
