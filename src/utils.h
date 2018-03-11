#ifndef BEACHMAT_UTILS_H
#define BEACHMAT_UTILS_H

#include "beachmat.h"

namespace beachmat { 

std::string make_to_string(const Rcpp::RObject&);

void throw_custom_error(const std::string&, const std::string&, const std::string&);

// Class checking.

std::string get_class(const Rcpp::RObject&);

Rcpp::RObject get_safe_slot(const Rcpp::RObject&, const std::string&);

std::string check_Matrix_class (const Rcpp::RObject&, const std::string&);

// Type checking.

std::string translate_type(int);

int find_sexp_type (const Rcpp::RObject&);

// Delayed Array conversion utilities.

Rcpp::RObject realize_delayed_array(const Rcpp::RObject&);

void check_DelayedMatrix(const Rcpp::RObject&);

bool only_delayed_coord_changes (const Rcpp::RObject&);

Rcpp::RObject extract_seed(const Rcpp::RObject&, const std::vector<std::string>&);

bool is_pristine_delayed_array(const Rcpp::RObject&);

// Matrix type enumeration.

enum matrix_type { SIMPLE, HDF5, SPARSE, RLE, PSYMM, DENSE, DELAYED };

}

#endif
