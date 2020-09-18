#include "beachmat3/beachmat.h"

template <class V, typename T = typename V::stored_type>
Rcpp::RObject get_sparse_row_slice0(Rcpp::RObject mat, Rcpp::IntegerVector order, 
    Rcpp::IntegerVector starts, Rcpp::IntegerVector ends) 
{
    auto ptr = beachmat::read_lin_sparse_block(mat);
    std::vector<int> work_i(ptr->get_ncol());
    std::vector<T> work_x(ptr->get_ncol());
    std::map<std::pair<int, int>, T> store;

    for (auto o : order) {
        int curstart = starts[o];
        int curend = ends[o];

        auto stuff = ptr->get_row(o, work_x.data(), work_i.data(), curstart, curend);
        for (size_t j = 0; j < stuff.n; ++j) {
            store[std::make_pair(o, stuff.i[j])] = stuff.x[j];
        }
    }

    return beachmat::as_gCMatrix<V>(ptr->get_nrow(), ptr->get_ncol(), store); 
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject get_sparse_row_slice(Rcpp::RObject mat, Rcpp::IntegerVector order, 
    Rcpp::IntegerVector starts, Rcpp::IntegerVector ends, int mode) 
{
    if (mode == 0) {
        return get_sparse_row_slice0<Rcpp::LogicalVector>(mat, order, starts, ends);
    } else {
        return get_sparse_row_slice0<Rcpp::NumericVector>(mat, order, starts, ends);
    }
}

template <class V, typename T = typename V::stored_type>
Rcpp::RObject get_sparse_row0(Rcpp::RObject mat, Rcpp::IntegerVector order) {
    auto ptr = beachmat::read_lin_sparse_block(mat);
    std::vector<int> work_i(ptr->get_ncol());
    std::vector<T> work_x(ptr->get_ncol());

    V store(ptr->get_nnzero());
    auto sIt = store.begin();
    Rcpp::IntegerVector newi(ptr->get_nnzero());
    auto iIt = newi.begin();

    for (auto o : order) {
        auto stuff = ptr->get_row(o, work_x.data(), work_i.data());
        for (size_t j = 0; j < stuff.n; ++j, ++sIt, ++iIt) {
            *iIt = stuff.i[j];
            *sIt = stuff.x[j];
        }
    }

    Rcpp::RObject output = beachmat::as_gCMatrix<V>(mat, store);
    output.slot("i") = newi;
    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject get_sparse_row(Rcpp::RObject mat, Rcpp::IntegerVector order, int mode) {
    if (mode == 0) {
        return get_sparse_row0<Rcpp::LogicalVector>(mat, order);
    } else {
        return get_sparse_row0<Rcpp::NumericVector>(mat, order);
    }
}
