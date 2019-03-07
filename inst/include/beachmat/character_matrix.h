#ifndef BEACHMAT_CHARACTER_MATRIX_H
#define BEACHMAT_CHARACTER_MATRIX_H

#include "Rcpp.h"

#include "simple_reader.h"
#include "dense_reader.h"
#include "delayed_reader.h"
#include "unknown_reader.h"
#include "external_reader.h"

#include "simple_writer.h"
#include "output_param.h"

#include "utils.h"

#include <memory>
#include <vector>

namespace beachmat { 

/*********
 * INPUT *
 *********/

std::unique_ptr<character_matrix> create_character_matrix_internal(const Rcpp::RObject&, bool); 

/* Virtual base class for character matrices. */

class character_matrix {
public:    
    character_matrix() = default;
    virtual ~character_matrix() = default;
    character_matrix(const character_matrix&) = default;
    character_matrix& operator=(const character_matrix&) = default;
    character_matrix(character_matrix&&) = default;
    character_matrix& operator=(character_matrix&&) = default;

    virtual size_t get_nrow() const=0;
    virtual size_t get_ncol() const=0;
    
    // Basic getters.
    void get_col(size_t c, Rcpp::StringVector::iterator out) { 
        get_col(c, out, 0, get_nrow());
        return;
    }
    virtual void get_col(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    void get_row(size_t r, Rcpp::StringVector::iterator out) { 
        get_row(r, out, 0, get_ncol());
        return;
    }
    virtual void get_row(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    virtual Rcpp::String get(size_t, size_t)=0;

    // Multi getters.
    void get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out) {
        get_cols(it, n, out, 0, get_nrow());
        return;
    }
    virtual void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    void get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out) {
        get_rows(it, n, out, 0, get_ncol());
        return;
    }
    virtual void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    // Specialized getters.
    Rcpp::StringVector::iterator get_const_col(size_t c, Rcpp::StringVector::iterator work) {
        return get_const_col(c, work, 0, get_nrow());
    }

    Rcpp::StringVector::iterator get_const_col(size_t c, Rcpp::StringVector::iterator work, size_t first, size_t last) {
        get_col(c, work, first, last);
        return work;
    }

    virtual Rcpp::StringVector::iterator get_const_col(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    const_col_indexed_info<Rcpp::StringVector> get_const_col_indexed(size_t c, Rcpp::StringVector::iterator work) {
        return get_const_col_indexed(c, work, 0, get_nrow());
    }

    const_col_indexed_info<Rcpp::StringVector> get_const_col_indexed(size_t c, Rcpp::StringVector::iterator work, size_t first, size_t last) {
        if (static_cast<size_t>(indices.size())!=this->get_nrow()) {
            indices=Rcpp::IntegerVector(this->get_nrow());
            std::iota(indices.begin(), indices.end(), 0); // populating with indices.
        }
        return const_col_indexed_info<Rcpp::StringVector>(last - first, indices.begin() + first, get_const_col(c, work, first, last));
    }

    virtual const_col_indexed_info<Rcpp::StringVector> get_const_col_indexed(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    // Other methods.
    virtual std::unique_ptr<character_matrix> clone() const=0;

    virtual Rcpp::RObject yield () const=0;
    virtual matrix_type get_matrix_type() const=0;

private:
    Rcpp::IntegerVector indices; // needed for get_const_col_indexed for non-sparse matrices.
};

/* Advanced character matrix template */

template<class RDR>
class general_character_matrix : public character_matrix {
public:    
    general_character_matrix(const Rcpp::RObject& incoming) : reader(incoming) {}
    ~general_character_matrix() = default;
    general_character_matrix(const general_character_matrix&) = default;
    general_character_matrix& operator=(const general_character_matrix&) = default;
    general_character_matrix(general_character_matrix&&) = default;
    general_character_matrix& operator=(general_character_matrix&&) = default;
  
    size_t get_nrow() const { return reader.get_nrow(); }
    size_t get_ncol() const { return reader.get_ncol(); }

    // Basic getters.
    void get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        reader.get_row(r, out, first, last);
        return;
    }
    void get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
        reader.get_col(c, out, first, last); 
        return;
    }

    Rcpp::String get(size_t r, size_t c) { return reader.get(r, c); }

    // Multi getters.
    void get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
        reader.get_rows(it, n, out, first, last);
        return;
    }
    void get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        reader.get_cols(it, n, out, first, last);
        return;
    }

    // Other methods.
    std::unique_ptr<character_matrix> clone() const { return std::unique_ptr<character_matrix>(new general_character_matrix(*this)); }

    Rcpp::RObject yield () const { return reader.yield(); }
    matrix_type get_matrix_type() const { return reader.get_matrix_type(); }
protected:
    RDR reader;
};

/* Simple character matrix */

using simple_character_precursor=general_character_matrix<simple_reader<Rcpp::String, Rcpp::StringVector> >;

class simple_character_matrix : public simple_character_precursor {
public:
    simple_character_matrix(const Rcpp::RObject& incoming) : simple_character_precursor (incoming) {}
    ~simple_character_matrix() = default;
    simple_character_matrix(const simple_character_matrix&) = default;
    simple_character_matrix& operator=(const simple_character_matrix&) = default;
    simple_character_matrix(simple_character_matrix&&) = default;
    simple_character_matrix& operator=(simple_character_matrix&&) = default;

    Rcpp::StringVector::iterator get_const_col(size_t c, Rcpp::StringVector::iterator work, size_t first, size_t last) {
        return reader.get_const_col(c, first, last);
    }

    std::unique_ptr<character_matrix> clone() const {
        return std::unique_ptr<character_matrix>(new simple_character_matrix(*this));
    }
};

/* DelayedMatrix */

typedef delayed_reader<Rcpp::String, Rcpp::StringVector, character_matrix> delayed_character_reader;

template<>
std::unique_ptr<character_matrix> delayed_character_reader::generate_seed(Rcpp::RObject incoming) {
    return create_character_matrix_internal(incoming, false);
}

using delayed_character_matrix=general_character_matrix<delayed_character_reader>;

/* Unknown matrix type */

using unknown_character_matrix=general_character_matrix<unknown_reader<Rcpp::String, Rcpp::StringVector> >;

/* External matrix type */

using external_character_precursor=general_character_matrix<external_reader<Rcpp::String, Rcpp::StringVector> >;

class external_character_matrix : public external_character_precursor {
public:
    external_character_matrix(const Rcpp::RObject& incoming) : external_character_precursor (incoming) {}
    ~external_character_matrix() = default;
    external_character_matrix(const external_character_matrix&) = default;
    external_character_matrix& operator=(const external_character_matrix&) = default;
    external_character_matrix(external_character_matrix&&) = default;
    external_character_matrix& operator=(external_character_matrix&&) = default;

    Rcpp::StringVector::iterator get_const_col(size_t c, Rcpp::StringVector::iterator work, size_t first, size_t last) {
        return reader.get_const_col(c, work, first, last);
    }

    const_col_indexed_info<Rcpp::StringVector> get_const_col_indexed(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        Rcpp::IntegerVector::iterator iIt;
        size_t nzero=this->reader.get_const_col_indexed(c, iIt, out, first, last);
        return const_col_indexed_info<Rcpp::StringVector>(nzero, iIt, out); 
    }

    std::unique_ptr<character_matrix> clone() const {
        return std::unique_ptr<character_matrix>(new external_character_matrix(*this));
    }
};

/* Dispatcher */

inline std::unique_ptr<character_matrix> create_character_matrix_internal(const Rcpp::RObject& incoming, bool delayed) { 
    if (incoming.isS4()) { 
        std::string ctype=get_class(incoming);
        if (ctype=="HDF5Matrix") {
            return std::unique_ptr<character_matrix>(new HDF5_character_matrix(incoming));
        } else if (delayed && ctype=="DelayedMatrix") { 
            return std::unique_ptr<character_matrix>(new delayed_character_matrix(incoming));
        } else if (has_external_support(incoming)) {
            return std::unique_ptr<character_matrix>(new external_character_matrix(incoming));
        }
        return std::unique_ptr<character_matrix>(new unknown_character_matrix(incoming));
    } 
    return std::unique_ptr<character_matrix>(new simple_character_matrix(incoming));
}

inline std::unique_ptr<character_matrix> create_character_matrix(const Rcpp::RObject& incoming) { 
    return create_character_matrix_internal(incoming, true);
}

/*********
 * INPUT *
 *********/

/* Virtual base class for character matrices. */

class character_output {
public:    
    character_output() = default;
    virtual ~character_output() = default;
    character_output(const character_output&) = default;
    character_output& operator=(const character_output&) = default;
    character_output(character_output&&) = default;
    character_output& operator=(character_output&&) = default;
    
    virtual size_t get_nrow() const=0;
    virtual size_t get_ncol() const=0;

    // Getters    
    void get_col(size_t c, Rcpp::StringVector::iterator out) { 
        get_col(c, out, 0, get_nrow());
        return;
    }
    virtual void get_col(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    void get_row(size_t r, Rcpp::StringVector::iterator out) { 
        get_row(r, out, 0, get_ncol());
        return;
    }
    virtual void get_row(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    virtual Rcpp::String get(size_t, size_t)=0;

    // Setters
    void set_col(size_t c, Rcpp::StringVector::iterator out) { 
        set_col(c, out, 0, get_nrow());
        return;
    }
    virtual void set_col(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    void set_row(size_t r, Rcpp::StringVector::iterator out) { 
        set_row(r, out, 0, get_ncol());
        return;
    }
    virtual void set_row(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    virtual void set(size_t, size_t, Rcpp::String)=0;

    virtual void set_col_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::StringVector::iterator)=0;

    virtual void set_row_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::StringVector::iterator)=0;

    // Other stuff.
    virtual Rcpp::RObject yield()=0;

    virtual std::unique_ptr<character_output> clone() const=0;

    virtual matrix_type get_matrix_type() const=0;
};

/* Simple character matrix */

class simple_character_output : public character_output {
public:
    simple_character_output(size_t nr, size_t nc) : writer(nr, nc) {}
    ~simple_character_output() = default;
    simple_character_output(const simple_character_output&) = default;
    simple_character_output& operator=(const simple_character_output&) = default;
    simple_character_output(simple_character_output&&) = default;
    simple_character_output& operator=(simple_character_output&&) = default;

    size_t get_nrow() const {
        return writer.get_nrow();
    }
    
    size_t get_ncol() const {
        return writer.get_ncol();
    }
    
    void get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
        writer.get_row(r, out, first, last);
        return;
    }
    
    void get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
        writer.get_col(c, out, first, last);
        return;
    }
    
    Rcpp::String get(size_t r, size_t c) {
        return writer.get(r, c);
    }
    
    void set_row(size_t r, Rcpp::StringVector::iterator in, size_t first, size_t last) { 
        writer.set_row(r, in, first, last);
        return;
    }
    
    void set_col(size_t c, Rcpp::StringVector::iterator in, size_t first, size_t last) { 
        writer.set_col(c, in, first, last);
        return;
    }
    
    void set(size_t r, size_t c, Rcpp::String in) {
        writer.set(r, c, in);
        return;
    }
    
    void set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
        writer.set_col_indexed(c, N, idx, val);
        return;
    }
    
    void set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
        writer.set_row_indexed(r, N, idx, val);
        return;
    }
    
    Rcpp::RObject yield() {
        return writer.yield();
    }
    
    std::unique_ptr<character_output> clone() const {
        return std::unique_ptr<character_output>(new simple_character_output(*this));
    }
    
    matrix_type get_matrix_type() const {
        return writer.get_matrix_type();
    }
private:
    simple_writer<Rcpp::String, Rcpp::StringVector> writer;
};

/* Dispatcher */

inline std::unique_ptr<character_output> create_character_output(int nrow, int ncol, const output_param& param) {
    switch (param.get_mode()) {
        case SIMPLE:
            return std::unique_ptr<character_output>(new simple_character_output(nrow, ncol));
        default:
            throw std::runtime_error("unsupported output mode for character matrices");
    }
}

}

#endif
