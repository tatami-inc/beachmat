#include "Rtatami.h"
#include "Rcpp.h"

#include "tatami/tatami.hpp"

#include <vector>
#include <cstddef>
#include <stdexcept>

#include "utils.h"

SEXP tatami_realize_svt_from_row(const tatami::NumericMatrix& mat, int threads) {
    const auto primary = mat.ncol();
    const auto secondary = mat.nrow();

    auto vbuffers = tatami::create_container_of_Index_size<std::vector<std::vector<double> > >(secondary);
    auto ibuffers = tatami::create_container_of_Index_size<std::vector<std::vector<int> > >(secondary);
    auto ranges = tatami::create_container_of_Index_size<std::vector<tatami::SparseRange<double, int> > >(secondary);

    tatami::parallelize([&](const int, const int start, const int len) -> void {
        auto ext = tatami::consecutive_extractor<true>(mat, true, start, len);
        auto vbuffer = tatami::create_container_of_Index_size<std::vector<double> >(primary);
        auto ibuffer = tatami::create_container_of_Index_size<std::vector<int> >(primary);
        for (auto i = start, end = start + len; i < end; ++i) {
            auto range = ext->fetch(vbuffer.data(), ibuffer.data());
            if (range.value == vbuffer.data()) {
                vbuffers[i].insert(vbuffers[i].end(), range.value, range.value + range.number);
                range.value = vbuffers[i].data();
            }
            if (range.index == ibuffer.data()) {
                ibuffers[i].insert(ibuffers[i].end(), range.index, range.index + range.number);
                range.index = ibuffers[i].data();
            }
            ranges[i] = range;
        }
    }, secondary, threads);

    auto primary_counts = tatami::create_container_of_Index_size<std::vector<I<decltype(secondary)> > >(primary);
    for (I<decltype(secondary)> s = 0; s < secondary; ++s) {
        const auto& sec_indices = ranges[s].index;
        const auto& sec_values = ranges[s].value;
        const auto num = ranges[s].number;
        for (I<decltype(num)> n = 0; n < num; ++n) {
            if (sec_values[n]) {
                primary_counts[sec_indices[n]] += 1; // addition must be safe as this cannot exceed dimension extents.
            }
        }
    }

    std::vector<Rcpp::IntegerVector> output_i;
    output_i.reserve(primary);
    std::vector<Rcpp::NumericVector> output_v;
    output_v.reserve(primary);
    tatami::cast_Index_to_container_size<Rcpp::IntegerVector>(secondary); // check that casts within the loop are safe.
    tatami::cast_Index_to_container_size<Rcpp::NumericVector>(secondary);
    for (I<decltype(primary)> p = 0; p < primary; ++p) {
        output_i.emplace_back(primary_counts[p]);
        output_v.emplace_back(primary_counts[p]);
    }

    std::fill(primary_counts.begin(), primary_counts.end(), 0);
    for (I<decltype(secondary)> s = 0; s < secondary; ++s) {
        const auto& sec_values = ranges[s].value;
        const auto& sec_indices = ranges[s].index;
        const auto num = ranges[s].number;
        for (I<decltype(num)> n = 0; n < num; ++n) {
            if (sec_values[n]) {
                const auto curp = sec_indices[n];
                auto& pos = primary_counts[curp];
                output_v[curp][pos] = sec_values[n];
                output_i[curp][pos] = s;
                ++pos;
            }
        }
    }

    Rcpp::List svt_out(primary);
    bool non_empty = false;
    for (I<decltype(primary)> p = 0; p < primary; ++p) {
        if (output_i[p].size() == 0) {
            svt_out[p] = R_NilValue;
            continue;
        }
        non_empty = true;

        bool all_ones = true; 
        for (auto v : output_v[p]) {
            if (v != 1) {
                all_ones = false;
                break;
            }
        }

        if (all_ones) {
            svt_out[p] = Rcpp::List::create(R_NilValue, std::move(output_i[p]));
        } else {
            svt_out[p] = Rcpp::List::create(std::move(output_v[p]), std::move(output_i[p]));
        }
    }

    Rcpp::S4 output("SVT_SparseMatrix");
    if (non_empty) {
        output.slot("SVT") = std::move(svt_out);
    } else {
        output.slot("SVT") = R_NilValue;
    }
    output.slot("type") = "double";
    output.slot("dim") = Rcpp::IntegerVector::create(secondary, primary);
    return output;
}

SEXP tatami_realize_svt_from_column(const tatami::NumericMatrix& mat, int threads) {
    const auto primary = mat.ncol();
    const auto secondary = mat.nrow();

    auto ext = tatami::consecutive_extractor<true>(mat, false, 0, primary);
    auto vbuffer = tatami::create_container_of_Index_size<std::vector<double> >(secondary);
    auto ibuffer = tatami::create_container_of_Index_size<std::vector<int> >(secondary);
    tatami::cast_Index_to_container_size<Rcpp::NumericVector>(secondary);
    tatami::cast_Index_to_container_size<Rcpp::IntegerVector>(secondary);

    Rcpp::List svt_out(primary);
    bool non_empty = false;
    for (I<decltype(primary)> p = 0; p < primary; ++p) {
        auto range = ext->fetch(vbuffer.data(), ibuffer.data());
        
        I<decltype(secondary)> num_non_zero = 0, num_ones = 0; 
        for (I<decltype(range.number)> n = 0; n < range.number; ++n) {
            num_non_zero += (range.value[n] != 0);
            num_ones += (range.value[n] == 1);
        }

        if (num_non_zero == 0) {
            svt_out[p] = R_NilValue;
            continue;
        }
        non_empty = true;

        Rcpp::NumericVector outv;
        Rcpp::IntegerVector outi;
        if (num_non_zero == range.number) {
            outi = Rcpp::IntegerVector(range.index, range.index + range.number);
            if (num_non_zero > num_ones) {
                outv = Rcpp::NumericVector(range.value, range.value + range.number);
            }
        } else {
            outi = Rcpp::IntegerVector(num_non_zero);
            if (num_non_zero == num_ones) {
                auto iIt = outi.begin();
                for (I<decltype(range.number)> n = 0; n < range.number; ++n) {
                    if (range.value[n]) {
                        *(iIt++) = range.index[n];
                    }
                }
            } else {
                outv = Rcpp::NumericVector(num_non_zero);
                auto iIt = outi.begin();
                auto vIt = outv.begin();
                for (I<decltype(range.number)> n = 0; n < range.number; ++n) {
                    if (range.value[n]) {
                        *(iIt++) = range.index[n];
                        *(vIt++) = range.value[n];
                    }
                }
            }
        }

        if (num_non_zero == num_ones) {
            svt_out[p] = Rcpp::List::create(R_NilValue, std::move(outi));
        } else {
            svt_out[p] = Rcpp::List::create(std::move(outv), std::move(outi));
        }
    }

    Rcpp::S4 output("SVT_SparseMatrix");
    if (non_empty) {
        output.slot("SVT") = std::move(svt_out);
    } else {
        output.slot("SVT") = R_NilValue;
    }
    output.slot("type") = "double";
    output.slot("dim") = Rcpp::IntegerVector::create(secondary, primary);
    return output;
}

SEXP tatami_realize_dgc(const tatami::NumericMatrix& mat, int threads) {
    tatami::RetrieveFragmentedSparseContentsOptions options;
    options.two_pass = false;
    options.num_threads = threads;
    auto frag = tatami::retrieve_fragmented_sparse_contents<double, int>(mat, false, options);
    const auto& store_v = frag.value;
    const auto& store_i = frag.index;

    const auto primary = mat.ncol();
    Rcpp::IntegerVector output_p(sanisizer::sum<I<decltype(std::declval<Rcpp::IntegerVector>().size())> >(primary, 1));
    for (I<decltype(primary)> p = 0; p < primary; ++p) {
        output_p[p + 1] = sanisizer::sum<int>(output_p[p], store_v[p].size());
    }

    const auto last_p = output_p[primary];
    auto output_v = tatami::create_container_of_Index_size<Rcpp::NumericVector>(last_p);
    auto output_i = tatami::create_container_of_Index_size<Rcpp::IntegerVector>(last_p);
    I<decltype(last_p)> offset = 0;
    for (I<decltype(primary)> p = 0; p < primary; ++p) {
        std::copy(store_v[p].begin(), store_v[p].end(), output_v.begin() + offset);
        std::copy(store_i[p].begin(), store_i[p].end(), output_i.begin() + offset);
        offset += store_v[p].size();
    }

    Rcpp::S4 output("dgCMatrix");
    output.slot("x") = std::move(output_v);
    output.slot("i") = std::move(output_i);
    output.slot("p") = std::move(output_p);
    output.slot("Dim") = Rcpp::IntegerVector::create(mat.nrow(), mat.ncol());
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP tatami_realize(SEXP raw_input, bool as_svt, int threads) {
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;

    if (!shared->sparse()) {
        auto output = create_output_matrix(shared->nrow(), shared->ncol());
        tatami::convert_to_dense(shared.get(), false, static_cast<double*>(output.begin()), threads);
        return output;
    } else if (as_svt) {
        if (shared->prefer_rows()) {
            return tatami_realize_svt_from_row(*shared, threads);
        } else {
            return tatami_realize_svt_from_column(*shared, threads);
        }
    } else {
        return tatami_realize_dgc(*shared, threads);
    }
}
