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
    Rcpp::List store(3);

    if (raw_p.sexp_type() != INTSXP) {
        throw std::runtime_error("'p' vector should be integer");
    }
    Rcpp::IntegerVector p(raw_p);
    store[0] = p;

    if (raw_i.sexp_type() != INTSXP) {
        throw std::runtime_error("'i' vector should be integer");
    }
    Rcpp::IntegerVector i(raw_i);
    store[1] = i;

    if (raw_x.sexp_type() == LGLSXP) {
        Rcpp::LogicalVector x(raw_x);
        store[2] = x;
        output->ptr.reset(store_sparse_matrix<int>(std::move(x), std::move(i), std::move(p), nrow, ncol, byrow));
    } else if (raw_x.sexp_type() == REALSXP) {
        Rcpp::NumericVector x(raw_x);
        store[2] = x;
        output->ptr.reset(store_sparse_matrix<double>(std::move(x), std::move(i), std::move(p), nrow, ncol, byrow));
    } else {
        throw std::runtime_error("'x' vector should be integer or real");
    }

    output->original = store; // holding references to all R objects created here, to avoid GC.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_SVT_SparseMatrix(int nr, int nc, Rcpp::RObject seed) {
    auto output = Rtatami::new_BoundNumericMatrix();

    std::vector<tatami::ArrayView<int> > indices(nc, tatami::ArrayView<int>(NULL, 0));
    std::vector<tatami::ArrayView<double> > values_d;
    std::vector<tatami::ArrayView<int> > values_i;

    std::string type = Rcpp::as<std::string>(seed.slot("type"));
    bool use_double = (type == "double");
    if (use_double) {
        values_d.resize(nc, tatami::ArrayView<double>(NULL, 0));
    } else if (type == "integer" || type == "logical") {
        values_i.resize(nc, tatami::ArrayView<int>(NULL, 0));
    } else {
        throw std::runtime_error("unsupported type '" + type + "' for a SVT_SparseMatrix");
    }

    // Setting up buffers of all-1 values so that we can create views on it
    // whenever we encounter lacunar leaf nodes.
    Rcpp::IntegerVector alloc_i;
    Rcpp::NumericVector alloc_d;

    // Storing everything that we take a view on to avoid GC. It may not be
    // enough to store a reference to the top-level seed, as there might be
    // some ALTREP magic happening under the hood to realize each vector.
    Rcpp::List store_i(nc), store_v(nc);

    tatami_r::parse_SVT_SparseMatrix(seed, [&](int c, const Rcpp::IntegerVector& curindices, bool all_ones, const auto& curvalues) {
        indices[c] = tatami::ArrayView<int>(static_cast<const int*>(curindices.begin()), curindices.size());
        store_i[c] = curindices;

        if (all_ones) {
            if (use_double) {
                if (static_cast<int>(alloc_d.size()) != nr) {
                    alloc_d = Rcpp::IntegerVector(nr);
                    std::fill(alloc_d.begin(), alloc_d.end(), 1);
                }
                values_d[c] = tatami::ArrayView<double>(static_cast<const double*>(alloc_d.begin()), curindices.size());
            } else {
                if (static_cast<int>(alloc_i.size()) != nr) {
                    alloc_i = Rcpp::IntegerVector(nr);
                    std::fill(alloc_i.begin(), alloc_i.end(), 1);
                }
                values_i[c] = tatami::ArrayView<int>(static_cast<const int*>(alloc_i.begin()), curindices.size());
            }
        } else {
            typedef tatami::ElementType<decltype(curvalues)> StoredValue;
            constexpr bool is_int = std::is_same<StoredValue, int>::value;
            if (is_int == use_double) {
                throw std::runtime_error("unexpected value vector type for a SVT_SparseMatrix of type '" + type + "'");
            }

            if constexpr(is_int) {
                values_i[c] = tatami::ArrayView<int>(static_cast<const int*>(curvalues.begin()), curvalues.size());
            } else {
                values_d[c] = tatami::ArrayView<double>(static_cast<const double*>(curvalues.begin()), curvalues.size());
            }
            store_v[c] = curvalues;
        }
    });

    if (use_double) {
        output->ptr.reset(new tatami::FragmentedSparseColumnMatrix<double, int, decltype(values_d), decltype(indices)>(nr, nc, std::move(values_d), std::move(indices)));
    } else {
        output->ptr.reset(new tatami::FragmentedSparseColumnMatrix<double, int, decltype(values_i), decltype(indices)>(nr, nc, std::move(values_i), std::move(indices)));
    }
    output->original = Rcpp::List::create(store_i, store_v, alloc_i, alloc_d);

    return output;
}
