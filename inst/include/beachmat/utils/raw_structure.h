#ifndef BEACHMAT_RAW_STRUCTURE_H
#define BEACHMAT_RAW_STRUCTURE_H

#include "Rcpp.h"
#include "copyable_vector.h"

namespace beachmat {

template<class V>
class raw_structure {
public:
    raw_structure(size_t nv=0, size_t ns=0) : values(nv), structure(ns), 
        values_start(values.vec.begin()), structure_start(structure.vec.begin()) {}

    ~raw_structure()=default;
    
    raw_structure(const raw_structure& in) : n(in.n), 
        use_own(in.use_own), structure(in.structure), values(in.values) 
    {
        initialize_iterators(in);
        return;
    }

    raw_structure& operator=(const raw_structure& in) {
        n=in.n; 
        use_own=in.use_own;
        structure=in.structure;
        values=in.values;
        initialize_iterators(in);
    }

    raw_structure(raw_structure&&)=default;
    raw_structure& operator=(raw_structure&&)=default;

    size_t get_n() const { return n; }
    Rcpp::IntegerVector::iterator get_structure_start() const { return structure_start; }
    typename V::iterator get_values_start() const { return values_start; }

    size_t& get_n() { return n; }
    Rcpp::IntegerVector::iterator& get_structure_start() { return structure_start; }
    typename V::iterator& get_values_start() { return values_start; }

    // Not for external use. If the vector ends up being reassigned 
    // (unlikely; you should prefer to use the constructor for choosing the size),
    // it is imporant that the iterators are re-initialized!
    Rcpp::IntegerVector& get_structure() { return structure.vec; }
    V& get_values() { return values.vec; }
private:
    bool use_own=true;
    copyable_holder<Rcpp::IntegerVector> structure;
    copyable_holder<V> values;

    size_t n=0;
    Rcpp::IntegerVector::iterator structure_start;
    typename V::iterator values_start;

    void initialize_iterators(const raw_structure& in) {
        if (use_own) {
            structure_start=structure.vec.begin();
            values_start=values.vec.begin();
        } else {
            structure_start=in.structure_start;
            values_start=in.values_start;
        }
    }
};

}

#endif
