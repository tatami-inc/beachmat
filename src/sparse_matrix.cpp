#include "Rtatami.h"

#include <cstdint>
#include <limits>
#include <type_traits>

template<typename T_, typename XVector_>
tatami::NumericMatrix* store_sparse_matrix(XVector_ x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol, bool byrow) {
    tatami::ArrayView<T_> x_view(static_cast<const T_*>(x.begin()), x.size());
    tatami::ArrayView<int> i_view(static_cast<const int*>(i.begin()), i.size());
    tatami::ArrayView<int> p_view(static_cast<const int*>(p.begin()), p.size());
    typedef tatami::CompressedSparseMatrix<double, int, decltype(x_view), decltype(i_view), decltype(p_view)> SparseMat;
    return new SparseMat(nrow, ncol, std::move(x_view), std::move(i_view), std::move(p_view), byrow, false);
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_sparse_matrix(Rcpp::RObject raw_x, Rcpp::RObject raw_i, Rcpp::RObject raw_p, int nrow, int ncol, bool byrow) {
    auto output = Rtatami::new_BoundNumericMatrix();
    output->original = Rcpp::List::create(raw_x, raw_i, raw_p); // holding references to all R objects used here.

    if (raw_p.sexp_type() != INTSXP) {
        throw std::runtime_error("'p' vector should be integer");
    }
    Rcpp::IntegerVector p(raw_p);

    if (raw_i.sexp_type() != INTSXP) {
        throw std::runtime_error("'i' vector should be integer");
    }
    Rcpp::IntegerVector i(raw_i);

    if (raw_x.sexp_type() == LGLSXP) {
        Rcpp::LogicalVector x(raw_x);
        output->ptr.reset(store_sparse_matrix<int>(std::move(x), std::move(i), std::move(p), nrow, ncol, byrow));
    } else if (raw_x.sexp_type() == REALSXP) {
        Rcpp::NumericVector x(raw_x);
        output->ptr.reset(store_sparse_matrix<double>(std::move(x), std::move(i), std::move(p), nrow, ncol, byrow));
    } else {
        throw std::runtime_error("'x' vector should be integer or real");
    }

    return output;
}

template<class XVector_, SEXPTYPE desired_sexp_>
tatami::NumericMatrix* parse_SVT_SparseMatrix_internal(int NR, int NC, Rcpp::RObject seed) {
    Rcpp::List svt = seed.slot("SVT");
    if (svt.size() != NC) {
        throw std::runtime_error("'SVT' slot in a SVT_SparseMatrix object should have length equal to the number of columns");
    }

    std::vector<tatami::ArrayView<int> > indices;
    typedef typename std::remove_reference<decltype(std::declval<XVector_>()[0])>::type StoredValue;
    std::vector<tatami::ArrayView<StoredValue> > values;
    indices.reserve(NC);
    values.reserve(NC);

    for (int c = 0; c < NC; ++c) {
        Rcpp::List inner = svt[c];
        if (inner.size() != 2) {
            throw std::runtime_error("each entry of the 'SVT' slot of a SVT_SparseMatrix object should be a list of length 2");
        }

        // Verify type to ensure that we're not making a view on a temporary array.
        Rcpp::RObject first = inner[0];
        if (first.sexp_type() != INTSXP) {
            throw std::runtime_error("first entry of each element of the 'SVT' slot in a SVT_SparseMatrix object should be an integer vector");
        }
        Rcpp::IntegerVector curindices(first);

        // Check for index contents is done inside the fragmented constructor.
        Rcpp::RObject second(inner[1]);
        if (second.sexp_type() != desired_sexp_) {
            throw std::runtime_error("second entry of an element of the 'SVT' slot in a SVT_SparseMatrix object has an unexpected type");
        }
        XVector_ curvalues(second);

        indices.emplace_back(static_cast<const int*>(curindices.begin()), curindices.size());
        values.emplace_back(static_cast<const StoredValue*>(curvalues.begin()), curvalues.size());
    }

    return new tatami::FragmentedSparseColumnMatrix<double, int, decltype(values), decltype(indices)>(NR, NC, std::move(values), std::move(indices));
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_SVT_SparseMatrix(int nr, int nc, Rcpp::RObject seed) {
    auto output = Rtatami::new_BoundNumericMatrix();
    output->original = seed;

    std::string type = Rcpp::as<std::string>(seed.slot("type"));
    if (type == "double") {
        output->ptr.reset(parse_SVT_SparseMatrix_internal<Rcpp::NumericVector, REALSXP>(nr, nc, seed));
    } else if (type == "integer") {
        output->ptr.reset(parse_SVT_SparseMatrix_internal<Rcpp::IntegerVector, INTSXP>(nr, nc, seed));
    } else if (type == "logical") {
        output->ptr.reset(parse_SVT_SparseMatrix_internal<Rcpp::LogicalVector, LGLSXP>(nr, nc, seed));
    } else {
        throw std::runtime_error("unsupported type '" + type + "' for a SVT_SparseMatrix");
    }

    return output;
}
