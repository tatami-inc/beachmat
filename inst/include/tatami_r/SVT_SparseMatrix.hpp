#ifndef TATAMI_R_SVT_SPARSE_MATRIX_HPP
#define TATAMI_R_SVT_SPARSE_MATRIX_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include <type_traits>

namespace tatami_r { 

template<typename Data_ = double, typename Index_ = int, class InputObject_, SEXPTYPE desired_sexp_>
Parsed<Data_, Index_> parse_SVT_SparseMatrix_internal(Rcpp::RObject seed) {
    auto dims = parse_dims(seed.slot("dim"));
    int NR = dims.first;
    int NC = dims.second;

    Rcpp::List svt = seed.slot("SVT");
    if (svt.size() != NC) {
        auto ctype = get_class_name(seed);
        throw std::runtime_error(std::string("'SVT' slot in a ") + ctype + " object should have length equal to the number of columns");
    }

    std::vector<tatami::ArrayView<int> > indices;
    typedef typename std::remove_reference<decltype(std::declval<InputObject_>()[0])>::type StoredValue;
    std::vector<tatami::ArrayView<StoredValue> > values;
    indices.reserve(NC);
    values.reserve(NC);

    for (int c = 0; c < NC; ++c) {
        Rcpp::List inner = svt[c];
        if (inner.size() != 2) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("each entry of the 'SVT' slot of a " + ctype + " object should be a list of length 2");
        }

        // Verify type to ensure that we're not making a view on a temporary array.
        Rcpp::RObject first = inner[0];
        if (first.sexp_type() != INTSXP) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("first entry of each element of the 'SVT' slot in a " + ctype + " object should be an integer vector");
        }
        Rcpp::IntegerVector curindices(first);

        // Check for index contents is done inside the fragmented constructor.
        Rcpp::RObject second(inner[1]);
        if (second.sexp_type() != desired_sexp_) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("second entry of an element of the 'SVT' slot in a " + ctype + " object has an unexpected type");
        }
        InputObject_ curvalues(second);

        indices.emplace_back(static_cast<const int*>(curindices.begin()), curindices.size());
        values.emplace_back(static_cast<const StoredValue*>(curvalues.begin()), curvalues.size());
    }

    Parsed<Data_, Index_> output;
    output.contents = seed;
    output.matrix.reset(new tatami::FragmentedSparseColumnMatrix<Data_, Index_, decltype(values), decltype(indices)>(NR, NC, std::move(values), std::move(indices)));
    return output;
}

template<typename Data_ = double, typename Index_ = int>
Parsed<Data_, Index_> parse_SVT_SparseMatrix(Rcpp::RObject seed) {
    std::string type = Rcpp::as<std::string>(seed.slot("type"));

    Parsed<Data_, Index_> output;
    if (type == "double") {
        output = parse_SVT_SparseMatrix_internal<Data_, Index_, Rcpp::NumericVector, REALSXP>(seed);
    } else if (type == "integer") {
        output = parse_SVT_SparseMatrix_internal<Data_, Index_, Rcpp::IntegerVector, INTSXP>(seed);
    } else if (type == "logical") {
        output = parse_SVT_SparseMatrix_internal<Data_, Index_, Rcpp::LogicalVector, LGLSXP>(seed);
    } else {
        auto ctype = get_class_name(seed);
        throw std::runtime_error("unsupported type '" + type + "' for a " + ctype + "object");
    }

    return output;
}

}

#endif
