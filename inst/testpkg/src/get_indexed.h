#ifndef BEACHTEST_GET_INDEXED_H
#define BEACHTEST_GET_INDEXED_H
#include "beachtest.h"
#include <algorithm>

Rcpp::IntegerVector spawn_indices(size_t n) {
    Rcpp::IntegerVector out(n);
    std::iota(out.begin(), out.end(), 0);
    return out;
}

template <class T, class O, class M>  // M is automatically deduced.
O get_indexed_all (M ptr, Rcpp::IntegerVector ordering) {
    const size_t& nrows=ptr->get_nrow();
    O output(nrows, ordering.size());

    const bool use_sparse=(ptr->col_raw_type()=="sparse");
    auto raws=ptr->set_up_raw();
    T alternative(nrows);
    size_t c=0;

    auto vals=alternative.begin();
    size_t N=nrows;
    Rcpp::IntegerVector counter(spawn_indices(nrows));
    auto idx=counter.begin();

    for (auto o : ordering) {
        if (use_sparse) {
            ptr->get_col_raw(o-1, raws);
            N=raws.get_n();
            idx=raws.get_structure_start();
            vals=raws.get_values_start();
        } else {
            ptr->get_col(o-1, vals);
        }

        auto outcol=output.column(c);
        for (size_t i=0; i<N; ++i) {
            outcol[*(idx+i)] = *(vals+i);
        }
        ++c;
    }
    return output;
}

template <class T, class O, class M>  
O get_indexed_slice (M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerVector rows) {
    if (rows.size()!=2) { 
        throw std::runtime_error("'rows' should be an integer vector of length 2"); 
    }
    const int rstart=rows[0]-1, rend=rows[1];
    const int out_nrows=rend-rstart;
    O output(out_nrows, ordering.size());

    const bool use_sparse=(ptr->col_raw_type()=="sparse");
    auto raws=ptr->set_up_raw();
    T alternative(out_nrows);
    size_t c=0;

    auto vals=alternative.begin();
    Rcpp::IntegerVector counter(spawn_indices(ptr->get_nrow()));

    for (auto o : ordering) {
        size_t N;
        Rcpp::IntegerVector::iterator idx;

        if (use_sparse) { 
            ptr->get_col_raw(o-1, raws, rstart, rend);
            N=raws.get_n();
            idx=raws.get_structure_start();
            vals=raws.get_values_start();
        } else {
            ptr->get_col(o-1, vals, rstart, rend);
            idx=counter.begin()+rstart;
            N=rend-rstart;
        }

        auto curcol=output.column(c);
        for (size_t i=0; i<N; ++i) {
            curcol[*(idx + i) - rstart] = *(vals+i);
        }
        ++c;
    }

    return output;
}

template <class T, class M>  
Rcpp::List get_indexed_varslice (M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerMatrix rows) {
    if (rows.ncol()!=2) { 
        throw std::runtime_error("'rows' should be an integer matrix with two columns"); 
    }
    if (rows.nrow()!=ordering.size()) {
        throw std::runtime_error("'nrow(rows)' should be equal to 'length(ordering)'");
    }
    Rcpp::List output(ordering.size());
    const size_t nrows=ptr->get_nrow();

    const bool use_sparse=(ptr->col_raw_type()=="sparse");
    auto raws=ptr->set_up_raw();
    T alternative(nrows);
    size_t c=0;

    auto vals=alternative.begin();
    Rcpp::IntegerVector counter(spawn_indices(nrows));

    for (auto o : ordering) {
        auto cur_bounds=rows.row(c);
        int left=cur_bounds[0]-1, right=cur_bounds[1];
        size_t N;
        Rcpp::IntegerVector::iterator idx;

        if (use_sparse) {
            ptr->get_col_raw(o-1, raws, left, right);
            N=raws.get_n();
            idx=raws.get_structure_start();
            vals=raws.get_values_start();
        } else {
            ptr->get_col(o-1, vals, left, right);
            idx=counter.begin()+left;
            N=right-left;
        }

        T out(right-left);
        for (size_t i=0; i<N; ++i) {
            out[*(idx + i) - left] = *(vals+i);
        }
        output[c]=out;
        ++c;
    }

    return output;
}

#endif
