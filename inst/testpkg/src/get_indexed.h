#ifndef BEACHTEST_GET_INDEXED_H
#define BEACHTEST_GET_INDEXED_H
#include "beachtest.h"

template <class T, class O, class M>  // M is automatically deduced.
O get_indexed_all (M ptr, Rcpp::IntegerVector ordering) {
    const size_t& nrows=ptr->get_nrow();
    O output(nrows, ordering.size());

    size_t c=0;
    for (auto o : ordering) {
        auto info=ptr->get_const_col_indexed(o-1);
        size_t N=std::get<0>(info);
        auto idx=std::get<1>(info);
        auto vals=std::get<2>(info);

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
    const int nrows=rend-rstart;    

    O output(nrows, ordering.size());
    size_t c=0;
    for (auto o : ordering) {
        auto info=ptr->get_const_col_indexed(o-1, rstart, rend);
        size_t N=std::get<0>(info);
        auto idx=std::get<1>(info);
        auto vals=std::get<2>(info);

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

    size_t c=0;
    for (auto o : ordering) {
        auto cur_bounds=rows.row(c);
        int left=cur_bounds[0]-1, right=cur_bounds[1];

        auto info=ptr->get_const_col_indexed(o-1, left, right);
        size_t N=std::get<0>(info);
        auto idx=std::get<1>(info);
        auto vals=std::get<2>(info);

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
