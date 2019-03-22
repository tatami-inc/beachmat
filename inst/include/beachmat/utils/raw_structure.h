#ifndef BEACHMAT_RAW_STRUCTURE_H
#define BEACHMAT_RAW_STRUCTURE_H

#include "Rcpp.h"
#include "copyable_vector.h"

namespace beachmat {

template<class V>
class raw_structure {
public:
    raw_structure(size_t nv=0, bool ownval=false, size_t ns=0, bool ownstruct=false) :
        values(nv, ownval), structure(ns, ownstruct) {}

    size_t get_n() const { return n; }
    Rcpp::IntegerVector::iterator get_structure_start() const { return structure.start; }
    typename V::iterator get_values_start() const { return values.start; }

    size_t& get_n() { return n; }
    Rcpp::IntegerVector::iterator& get_structure_start() { return structure.start; }
    typename V::iterator& get_values_start() { return values.start; }

    // Not for external use. If the vector ends up being reassigned 
    // (unlikely; you should prefer to use the constructor for choosing the size),
    // it is imporant that the iterators are re-initialized!
    Rcpp::IntegerVector& get_structure() { return structure.holder.vec; }
    V& get_values() { return values.holder.vec; }
private:
    size_t n=0;

    // Internal class for easier copying and moving.
    template<typename T>
    struct internal_vector {
        copyable_holder<T> holder;
        typename T::iterator start;
        bool use_own;

        internal_vector(size_t n, bool u) : holder(n), start(holder.vec.begin()), use_own(u) {}
        ~internal_vector()=default;

        // Either copying the iterator directly, if we don't use our own structure;
        // or resetting it to an appropriate position in our own structure.
        internal_vector(const internal_vector& in) : holder(in.holder), use_own(in.use_own) {
            initialize_iterator(in);
            return;
        }
        internal_vector& operator=(const internal_vector& in) {
            holder=in.holder;
            use_own=in.use_own;
            initialize_iterator(in);
            return *this;
        }
        void initialize_iterator(const internal_vector& in) {
            if (use_own) {
                start=holder.vec.begin(); //+ (in.start - in.holder.vec.begin());
            } else {
                start=in.start;
            }
            return;
        }

        // Copying the iterator in the move operator; 
        // do not use the automatically generated move, this is not correct.
        internal_vector(internal_vector&& in) : holder(std::move(in.holder)), use_own(std::move(in.use_own)) {
            initialize_iterator(in);
        }
        internal_vector& operator=(internal_vector&& in) {
            holder=std::move(in.holder);
            use_own=std::move(in.use_own);
            initialize_iterator(in);
            return *this;
        } 
    }; 

    internal_vector<V> values;
    internal_vector<Rcpp::IntegerVector> structure;
};

}

#endif
