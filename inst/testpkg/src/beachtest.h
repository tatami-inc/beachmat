#ifndef BEACHTEST_H
#define BEACHTEST_H

#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/logical_matrix.h"
#include "beachmat/character_matrix.h"

Rcpp::StringVector translate_class(beachmat::matrix_type);

#endif
