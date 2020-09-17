#ifndef BEACHMAT_CSPARSE_MATRIX_H
#define BEACHMAT_CSPARSE_MATRIX_H

/**
 * @file Csparse_reader.h
 *
 * Internal utilities and class definitions for processing compressed sparse column matrix representations.
 */

#include "Rcpp.h"

#include "dim_checker.h"
#include "utils.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace beachmat {

/**
 * Sparse index container, holding the number of non-zero elements extracted from a single row/column 
 * along with pointers to their indices and values.
 * For extracted rows, the indices refer to column positions, and vice versa for extracted columns.
 * 
 * @tparam TIT The type of the (`const`) random-access iterator pointing to the data values.
 * @tparam I The integer type of the index.
 */
template <typename TIT, typename I>
struct sparse_index {
    /**
     * Constructor for the `sparse_index`, setting its data members directly to the supplied values without too much fuss.
     *
     * @param _n See `n`.
     * @param _x See `x`.
     * @param _i See `i`.
     */
    sparse_index(size_t _n, TIT _x, const I* _i) : n(_n), x(_x), i(_i) {}

    /**
     * Number of non-zero elements.
     */
    size_t n; 

    /**
     * Iterator pointing to the set of non-zero values.
     * This should be incrementable up to `n`.
     */
    TIT x;

    /**  
     * Pointer to the set of indices of the non-zero values.
     * This should be incrementable up to `n`.
     */
    const I* i;
};

/**
 * @internal
 *
 * Transplant indices and values into their respective workspaces 
 * and use them to construct a new `sparse_index` object.
 * This is necessary for type conversions between the native value type and the expected value type.
 *
 * @tparam OUT Iterator to be stored in the new `sparse_index`.
 * This should be the `const`-type counterpart to `ALT`.
 * @tparam TIT `const` iterator of native type, i.e., equal to the original data values in the matrix.
 * @tparam ALT Non-`const` iterator to the workspace of the desired type.
 * @tparam I Integer type of the index.
 *
 * @param ref `sparse_index` containing pointers from a column (usually) of a sparse representation.
 * @param work_x Pointer to the workspace for the data values, usually of a different type to that of the native representation.
 * This should have space for at least `ref.n` values.
 * @param work_i Pointeger to the workspace for the indices of the non-zero values.
 * This should have space for at least `ref.n` values.
 * 
 * @return A `sparse_index` containing iterators for `work_x` and `work_i`.
 */
template<typename OUT, typename TIT, typename ALT, typename I>
inline sparse_index<OUT, int> transplant(sparse_index<TIT, I> ref, ALT work_x, I* work_i) {
    std::copy(ref.x, ref.x + ref.n, work_x);
    std::copy(ref.i, ref.i + ref.n, work_i);
    return sparse_index<OUT, int>(ref.n, work_x, work_i);
}

/**
 * @internal
 *
 * Core handler for compressed sparse column (CSC) matrices,
 * controlling extraction of elements from the columns and rows for use by other classes.
 *
 * @tparam TIT The type of the (`const`) random-access iterator pointing to the data values.
 * @tparam I The integer type of the index.
 * @tparam P The integer type of the column pointers.
 */
template <typename TIT, typename I, typename P>
class Csparse_core {
public:
    /** 
     * Trivial constructor.
     */
    Csparse_core() {};

    /**
     * Constructor where arguments are copied directly into their corresponding data members.
     *
     * @param _n Number of non-zero elements.
     * @param _x Iterator to the non-zero data values in the matrix.
     * This should have at least `n` addressable elements.
     * @param _i Pointer to the row indices of the non-zero data values in the matrix.
     * This should have at least `n` addressable elements.
     * @param _nr Number of rows in the CSC matrix.
     * @param _nc Number of columns in the CSC matrix.
     * @param _p Pointer to the array of column pointers, 
     * specifying the starting index of each column in `_x` and `_i`.
     * This should have at least `nc + 1` addressable elements.
     */
    Csparse_core(const size_t _n, TIT _x, const I* _i, const size_t _nr, const size_t _nc, const P* _p) : 
        n(_n), nr(_nr), nc(_nc), x(_x), i(_i), p(_p), currow(0), curstart(0), curend(nc) {}
   
    /**
     * Get all non-zero elements from a column. 
     *
     * @param c The index of the column to extract.
     * @param first The index of the first row of interest.
     * @param last The index of the first row that is _not_ of interest.
     *
     * @return A `sparse_index` containing pointers to the first non-zero element in `c` with row index no less than `first`.
     * The number of non-zero elements is set to all those in `[first, last)`. 
     */
    sparse_index<TIT, I> get_col(size_t c, size_t first, size_t last) {
        const auto pstart=p[c]; 
        auto iIt = i + pstart, 
             eIt = i + p[c+1]; 
        auto xIt = x + pstart;

        if (first) { // Jumping ahead if non-zero.
            auto new_iIt=std::lower_bound(iIt, eIt, first);
            xIt+=(new_iIt-iIt);
            iIt=new_iIt;
        } 

        if (last!=nr) { // Jumping to last element.
            eIt=std::lower_bound(iIt, eIt, last);
        }

        return sparse_index<TIT, I>(eIt - iIt, xIt, iIt); 
    }

    /**
     * The type of the values pointed to by `TIT`. 
     */
    typedef decltype(*std::declval<TIT>()) T;

    /**
     * Get all elements from a column, explicitly filling in the zeroes.
     *
     * @tparam ALT Iterator class for the workspace.
     *
     * @param c The index of the column to extract.
     * @param work A pointer or iterator to the workspace in which the column values are to be stored.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first row of interest.
     * @param last The index of the first row that is _not_ of interest.
     * @param empty Value corresponding to zero, almost always `0`.
     *
     * @return `work` is filled in with the contents of column `c` from rows `[first, last)`.
     * If no non-zero element exists, the corresponding entry of `work` is set to `empty`.
     */
    template <typename ALT = TIT>
    void get_col(size_t c, ALT work, size_t first, size_t last, T empty) {
        auto out = this->get_col(c, first, last);
        std::fill(work, work + last - first, empty);
        for (size_t v = 0; v < out.n; ++v, ++out.i, ++out.x) {
            *(work + *out.i - first) = *out.x;
        }
        return;       
    }

    /**
     * Get all elements from a row, explicitly filling in the zeroes.
     *
     * @tparam ALT Iterator class for the workspace.
     *
     * @param r The index of the row to extract.
     * @param work A pointer or iterator to the workspace in which the row values are to be stored.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of the first column that is _not_ of interest.
     * @param empty Value corresponding to zero, almost always `0`.
     *
     * @return `work` is filled in with the contents of row `r` from columns `[first, last)`.
     * If no non-zero element exists, the corresponding entry of `work` is set to `empty`.
     */
    template <typename ALT = TIT>
    void get_row(size_t r, ALT work, size_t first, size_t last, T empty) {
        update_indices(r, first, last);
        std::fill(work, work + last - first, empty);

        auto pIt = p + first + 1; // Points to first-past-the-end for each 'c'.
        for (size_t c = first; c < last; ++c, ++pIt, ++work) { 
            const int idex = indices[c];
            if (static_cast<P>(idex) != *pIt && static_cast<size_t>(i[idex]) == r) { 
                (*work) = *(x + idex); 
            }
        } 
        return;  
    }

    /**
     * Get all non-zero elements from a row.
     *
     * @tparam OUT Iterator class for the data values in the output `sparse_index`,
     * expected to correspond to a `const`-type counterpart to `ALT`.
     * @tparam ALT Iterator class for the workspace.
     *
     * @param r The index of the row to extract.
     * @param work_x A pointer or iterator to the workspace in which the non-zero row values are to be stored.
     * This should have at least `last - first` addressable elements.
     * @param work_i A pointer or iterator to the workspace in which the non-zero column indices are to be stored.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of the first column that is _not_ of interest.
     * @param empty Value corresponding to zero, almost usually `0`.
     *
     * @return A `sparse_index` containing pointers to the workspaces.
     * The number of non-zero elements is set to all those in `[first, last)`. 
     */
    template <typename OUT, typename ALT = TIT>
    sparse_index<OUT, I> get_row(size_t r, ALT work_x, I* work_i, size_t first, size_t last) {
        update_indices(r, first, last);

        auto pIt = p + first + 1; // Points to first-past-the-end for each 'c'.
        size_t counter = 0;

        for (size_t c = first; c < last; ++c, ++pIt) { 
            const int& idex = indices[c];
            if (idex != *pIt && static_cast<size_t>(i[idex]) == r) { 
                work_i[counter] = i[idex];
                *(work_x + counter) = *(x + idex);
                ++counter;
            }
        }

        return sparse_index<OUT, I>(counter, work_x, work_i);
    }
private:
    size_t n, nr, nc;

    TIT x;
    const I* i;
    const P* p;

    size_t currow, curstart, curend;
    std::vector<P> indices; 

    /**
     * Update the index to the last requested non-zero element in each column.
     * This is used to accelerate consecutive row queries by avoiding the need for a new binary search.
     * After one row extraction, we hold the index of the lower-bounded non-zero element for each column;
     * on the next request for a row, we check whether it can be satisfied by the next non-zero element in that column.
     *
     * @param r The requested row.
     * @param first The first column of interest.
     * @param last The first column that is not of interest.
     * 
     * @return `indices` is updated.
     * 
     */
    void update_indices(size_t r, size_t first, size_t last) {
        /* Initializing the indices upon the first request, assuming currow=0 based on initialization above.
         * This avoids using up space for the indices if we never do row access.
         */
        if (indices.size() != nr) {
            indices = std::vector<P>(p, p + nc);
        }

        /* If left/right slice are not equal to what is stored, we reset the indices,
         * so that the code below will know to recompute them. It's too much effort
         * to try to figure out exactly which columns need recomputing; just do them all.
         */
        if (first != curstart || last != curend) {
            std::copy(p, p + nc, indices.begin());
            currow=0;
        }

        /* entry of 'indices' for each column should contain the index of the first
         * element with row number not less than 'r'. If no such element exists, it
         * will contain the index of the first element of the next column.
         */
        if (r == currow) { 
            return; 
        } 

        const P* pIt = p + first;
        if (r == currow+1) {
            ++pIt; // points to the first-past-the-end element, at any given 'c'.
            for (size_t c=first; c<last; ++c, ++pIt) {
                P& curdex = indices[c];
                if (curdex != *pIt && static_cast<size_t>(i[curdex]) < r) { 
                    ++curdex;
                }
            }
        } else if (r+1 == currow) {
            for (size_t c=first; c<last; ++c, ++pIt) {
                P& curdex = indices[c];
                if (curdex != *pIt && static_cast<size_t>(i[curdex-1]) >= r) { 
                    --curdex;
                }
            }

        } else { 
            if (r > currow) {
                ++pIt; // points to the first-past-the-end element, at any given 'c'.
                for (size_t c = first; c < last; ++c, ++pIt) { 
                    indices[c] = std::lower_bound(i + indices[c], i + *pIt, r) - i;
                }
            } else { 
                for (size_t c = first; c < last; ++c, ++pIt) {
                    indices[c] = std::lower_bound(i + *pIt, i + indices[c], r) - i;
                }
            }
        }

        currow=r;
        return;
    }
};

/**
 * @internal
 *
 * Reader for `dgCMatrix` and `lgCMatrix` R objects.
 *
 * @tparam V The type of the `Rcpp::Vector` containing the data values.
 * @tparam TIT The type of the (`const`) random-access iterator pointing to the data values.
 */
template <class V, typename TIT = typename V::iterator>
class gCMatrix_reader : public dim_checker {
public:
    ~gCMatrix_reader() = default;
    gCMatrix_reader(const gCMatrix_reader&) = default;
    gCMatrix_reader& operator=(const gCMatrix_reader&) = default;
    gCMatrix_reader(gCMatrix_reader&&) = default;
    gCMatrix_reader& operator=(gCMatrix_reader&&) = default;

    /**
     * Constructor from an R object containing a `*gCMatrix` instance.
     * This implements a series of checks for the validity of the slots for the compressed column-sparse format.
     *
     * @param mat An R object containing a `*gCMatrix` instance.
     */
    gCMatrix_reader(Rcpp::RObject mat) { 
        this->fill_dims(get_safe_slot(mat, "Dim"));
        const size_t& NC=this->ncol;
        const size_t& NR=this->nrow;

        Rcpp::RObject temp_i=get_safe_slot(mat, "i");
        if (temp_i.sexp_type()!=INTSXP) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("'i' slot in a ") + ctype + " object should be integer"); 
        }
        i=temp_i;

        Rcpp::RObject temp_p=get_safe_slot(mat, "p");
        if (temp_p.sexp_type()!=INTSXP) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should be integer");
        }
        p=temp_p;

        Rcpp::RObject temp_x=get_safe_slot(mat, "x");
        if (temp_x.sexp_type()!=x.sexp_type()) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("'x' slot in a ") + ctype + " object should be " + translate_type(x.sexp_type()));
        }
        x=temp_x;

        if (x.size()!=i.size()) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("'x' and 'i' slots in a ") + ctype + " object should have the same length"); 
        }
        if (NC+1!=static_cast<size_t>(p.size())) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("length of 'p' slot in a ") + ctype + " object should be equal to 'ncol+1'"); 
        }
        if (p[0]!=0) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("first element of 'p' in a ") + ctype + " object should be 0"); 
        }
        if (p[NC]!=x.size()) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("last element of 'p' in a ") + ctype + " object should be 'length(x)'"); 
        }

        // Checking all the indices.
        auto pIt=p.begin();
        for (size_t px=0; px<NC; ++px) {
            const int& current=*pIt;
            if (current < 0) { 
                auto ctype = get_class_name(mat);
                throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should contain non-negative values"); 
            }
            if (current > *(++pIt)) { 
                auto ctype = get_class_name(mat);
                throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should be sorted"); 
            }
        }

        pIt=p.begin();
        for (size_t px=0; px<NC; ++px) {
            int left=*pIt; // Integers as that's R's storage type. 
            int right=*(++pIt)-1; // Not checking the last element, as this is the start of the next column.
            auto iIt=i.begin()+left;

            for (int ix=left; ix<right; ++ix) {
                const int& current=*iIt;
                if (current > *(++iIt)) {
                    auto ctype = get_class_name(mat);
                    throw std::runtime_error(std::string("'i' in each column of a ") + ctype + " object should be sorted");
                }
            }
        }

        for (auto iIt=i.begin(); iIt!=i.end(); ++iIt) {
            const int& curi=*iIt;
            if (curi<0 || static_cast<size_t>(curi)>=NR) {
                auto ctype = get_class_name(mat);
                throw std::runtime_error(std::string("'i' slot in a ") + ctype + " object should contain elements in [0, nrow)");
            }
        }

        core=Csparse_core<TIT, int, int>(i.size(), x.begin(), i.begin(), NR, NC, p.begin());
        return;                
    }

    /**
     * @copydoc Csparse_core::get_col(size_t, size_t, size_t)
     */
    sparse_index<TIT, int> get_col(size_t c, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        return core.get_col(c, first, last);
    }

    /**
     * @copydoc Csparse_core::get_row(size_t, ALT, I*, size_t, size_t)
     */
    template <typename OUT, typename ALT = TIT>
    sparse_index<OUT, int> get_row(size_t r, ALT work_x, int* work_i, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        return core.template get_row<OUT>(r, work_x, work_i, first, last);
    }

    /**
     * @copydoc Csparse_core::get_row(size_t, ALT, size_t, size_t, T)
     */
    template <typename ALT = TIT>
    ALT get_row(size_t r, ALT work, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        core.get_row(r, work, first, last, 0);
        return work;
    }

    /**
     * @copydoc Csparse_core::get_col(size_t, ALT, size_t, size_t, T)
     */
    template <typename ALT = TIT>
    ALT get_col(size_t c, ALT work, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        core.get_col(c, work, first, last, 0);
        return work;
    }

    /**
     * Get the number of non-zero elements in the object.
     */
    size_t get_nnzero () const { return x.size(); }

private:
    Rcpp::IntegerVector i, p;
    V x;
    Csparse_core<TIT, int, int> core;
};

/**
 * @internal
 *
 * Reader for `SparseArraySeed` R objects.
 *
 * @tparam V The type of the `Rcpp::Vector` containing the data values.
 * @tparam TIT The type of the (`const`) random-access iterator pointing to the data values.
 */
template <class V, typename TIT = typename V::iterator>
class SparseArraySeed_reader : public dim_checker {
private:
    struct sparse_triplet {
        sparse_triplet(int _i, int _j, size_t _ptr) : i(_i), j(_j), ptr(_ptr) {}
        int i;
        int j;
        size_t ptr;
    };
public:
    ~SparseArraySeed_reader() = default;
    SparseArraySeed_reader(const SparseArraySeed_reader&) = default;
    SparseArraySeed_reader& operator=(const SparseArraySeed_reader&) = default;
    SparseArraySeed_reader(SparseArraySeed_reader&&) = default;
    SparseArraySeed_reader& operator=(SparseArraySeed_reader&&) = default;

    /**
     * Constructor from an R object containing a `SparseArraySeed` instance.
     * This implements a series of checks for the consistency of the slots.
     * It will copy the row and column indices to make them zero-based,
     * and if necessary, sort them to create a compressed sparse column format.
     *
     * @param mat An R object containing a `SparseArraySeed` instance.
     */
    SparseArraySeed_reader(Rcpp::RObject seed) {
        this->fill_dims(get_safe_slot(seed, "dim"));
        const size_t& NC=this->ncol;
        const size_t& NR=this->nrow;
        p.resize(NC + 1);
    
        Rcpp::RObject temp_x=get_safe_slot(seed, "nzdata");
        if (temp_x.sexp_type()!=x.sexp_type()) { 
            auto ctype = get_class_name(seed);
            throw std::runtime_error(std::string("'nzdata' slot in a ") + ctype + " object should be " 
                + translate_type(x.sexp_type()));
        }
        x=temp_x;

        Rcpp::RObject temp_i=get_safe_slot(seed, "nzindex");
        if (temp_i.sexp_type() != INTSXP) { 
            auto ctype = get_class_name(seed);
            throw std::runtime_error(std::string("'nzindex' slot in a ") + ctype + " object should be integer"); 
        }
        Rcpp::IntegerMatrix temp_i2(temp_i);
        if (temp_i2.ncol() != 2) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error(std::string("'nzindex' slot in a ") + ctype + " object should have two columns"); 
        }
        const size_t nnz = temp_i2.nrow();
        if (nnz != x.size()) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error(std::string("incompatible 'nzindex' and 'nzdata' lengths in a ") + ctype + " object"); 
        }

        if (!nnz) {
            auto row_indices=temp_i2.column(0);
            auto col_indices=temp_i2.column(1);
            bool okay=true;
            {
                auto rowIt=row_indices.begin();
                auto colIt=col_indices.begin();

                for (size_t v = 0; v < nnz; ++v) {
                    auto lastR=*rowIt;
                    auto lastC=*colIt;
                    if (lastR <= 0 || lastR > NR || lastC <= 0 || lastC > NC) {
                        auto ctype = get_class_name(seed);
                        throw std::runtime_error(std::string("'nzindex' out of bounds in a ") + ctype + " object");
                    }

                    if (okay && v < nnz - 1) {
                        auto nextR=*(++rowIt);
                        auto nextC=*(++colIt);
                        if (lastC > nextC || (lastC == nextC && lastR > nextR)) {
                            okay = false;
                        }
                    }
                }
            }

            if (okay) {
                i = Rcpp::IntegerVector(row_indices.begin(), row_indices.end());
                for (auto& subi : i) { --subi; }

                auto colIt = col_indices.begin();
                auto colStart = colIt;
                for (int c = 1; c <= static_cast<int>(NC); ++c) {
                    // Technically it should be *colIt <= c+1 to get to 1-based
                    // indices, but this cancels out with a -1 because we want
                    // everything up to the _last_ column.
                    while (colIt != col_indices.end() && *colIt <= c) { 
                        ++colIt;
                    }
                    p[c] = colIt - colStart;
                }
            } else {
                std::vector<sparse_triplet> everything;
                everything.reserve(nnz);

                auto rowIt=row_indices.begin();
                auto colIt=col_indices.begin();
                auto xIt=x.begin();

                for (int v = 0; v < nnz; ++v, ++rowIt, ++colIt, ++xIt) {
                    everything.push_back(sparse_triplet(*rowIt, *colIt, v));
                }

                std::sort(everything.begin(), everything.end(), 
                    [] (const sparse_triplet& left, const sparse_triplet& right) -> bool {
                    if (left.j < right.j) {
                        return true;
                    } else if (left.j > right.j) {
                        return false;
                    } else if (left.i < right.i) {
                        return true;
                    } else if (left.i > right.i) {
                        return false;
                    } else if (left.ptr < right.ptr) {
                        return true;
                    } else {
                        return false;
                    }
                });

                i = Rcpp::IntegerVector(nnz);
                V clone_x(nnz);
                for (size_t v = 0; v < nnz; ++v) {
                    i[v] = everything[v].i - 1;
                    clone_x[v] = x[everything[v].ptr];
                }
                x = clone_x;

                auto eIt = everything.begin();
                auto eStart = eIt;
                for (int c = 1; c <= static_cast<int>(NC); ++c) {
                    while (eIt != everything.end() && eIt->j <= c) { // see above comments.
                        ++eIt;
                    }
                    p[c] = eIt - eStart;
                }
            }
        }

        core=Csparse_core<TIT, int, size_t>(nnz, x.begin(), i.begin(), NR, NC, p.data());
        return;
    }

    /**
     * @copydoc Csparse_core::get_col(size_t, size_t, size_t)
     */
    sparse_index<TIT, int> get_col(size_t c, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        return core.get_col(c, first, last);
    }

    /**
     * @copydoc Csparse_core::get_row(size_t, ALT, I*, size_t, size_t)
     */
    template <typename OUT, typename ALT = TIT>
    sparse_index<OUT, int> get_row(size_t r, ALT work_x, int* work_i, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        return core.template get_row<OUT>(r, work_x, work_i, first, last);
    }

    /**
     * @copydoc Csparse_core::get_row(size_t, ALT, size_t, size_t, T)
     */
    template <typename ALT = TIT>
    ALT get_row(size_t r, ALT work, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        core.get_row(r, work, first, last, 0);
        return work;
    }

    /**
     * @copydoc Csparse_core::get_col(size_t, ALT, size_t, size_t, T)
     */
    template <typename ALT = TIT>
    ALT get_col(size_t c, ALT work, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        core.get_col(c, work, first, last, 0);
        return work;
    }

    /**
     * Get the number of non-zero elements in the object.
     */
    size_t get_nnzero () const { return x.size(); }

private:
    Rcpp::IntegerVector i;
    V x;
    std::vector<size_t> p;
    Csparse_core<TIT, int, size_t> core;
};

}

#endif
