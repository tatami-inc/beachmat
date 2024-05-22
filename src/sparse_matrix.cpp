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

//[[Rcpp::export(rng=false)]]
SEXP initialize_SVT_SparseMatrix(int nr, int nc, Rcpp::RObject seed) {
    auto output = Rtatami::new_BoundNumericMatrix();

    std::vector<tatami::ArrayView<int> > indices;
    indices.reserve(nc);
    std::vector<tatami::ArrayView<double> > values_d;
    std::vector<tatami::ArrayView<int> > values_i;

    // Creating buffers of all-1 values so that we can create views on it,
    // whenever we encounter lacunar leaf nodes.
    Rcpp::IntegerVector alloc_i;
    Rcpp::NumericVector alloc_d;

    std::string type = Rcpp::as<std::string>(seed.slot("type"));
    bool use_double = false;
    if (type == "double") {
        use_double = true;
        values_d.reserve(nc);
    } else if (type == "integer" || type == "logical") {
        values_i.reserve(nc);
    } else {
        throw std::runtime_error("unsupported type '" + type + "' for a SVT_SparseMatrix");
    }

    tatami_r::parse_SVT_SparseMatrix(seed, [&](int c, const Rcpp::IntegerVector& curindices, bool all_ones, const auto& curvalues) {
        indices.emplace_back(static_cast<const int*>(curindices.begin()), curindices.size());

        typedef tatami::ElementType<decltype(curvalues)> StoredValue;
        constexpr bool is_int = std::is_same<StoredValue, int>::value;

        if (all_ones) {
            if constexpr(is_int) {
                if (static_cast<int>(alloc_i.size()) != nc) {
                    alloc_i = Rcpp::IntegerVector(nc);
                    std::fill(alloc_i.begin(), alloc_i.end(), 1);
                }
                values_i.emplace_back(static_cast<const int*>(alloc_i.begin()), curindices.size());

            } else {
                if (static_cast<int>(alloc_d.size()) != nc) {
                    alloc_d = Rcpp::IntegerVector(nc);
                    std::fill(alloc_d.begin(), alloc_d.end(), 1);
                }
                values_d.emplace_back(static_cast<const double*>(alloc_d.begin()), curindices.size());
            }

        } else {
            if (is_int == use_double) {
                throw std::runtime_error("unexpected value vector type for a SVT_SparseMatrix of type '" + type + "'");
            }
            if constexpr(is_int) {
                values_i.emplace_back(static_cast<const int*>(curvalues.begin()), curvalues.size());
            } else {
                values_d.emplace_back(static_cast<const double*>(curvalues.begin()), curvalues.size());
            }
        }
    });

    output->original = Rcpp::List::create(seed, alloc_i, alloc_d);
    if (use_double) {
        output->ptr.reset(new tatami::FragmentedSparseColumnMatrix<double, int, decltype(values_d), decltype(indices)>(nr, nc, std::move(values_d), std::move(indices)));
    } else {
        output->ptr.reset(new tatami::FragmentedSparseColumnMatrix<double, int, decltype(values_i), decltype(indices)>(nr, nc, std::move(values_i), std::move(indices)));
    }

    return output;
}
