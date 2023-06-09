// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// apply_delayed_binary_operation
SEXP apply_delayed_binary_operation(SEXP left_input, SEXP right_input, std::string op);
RcppExport SEXP _beachmat_apply_delayed_binary_operation(SEXP left_inputSEXP, SEXP right_inputSEXP, SEXP opSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type left_input(left_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type right_input(right_inputSEXP);
    Rcpp::traits::input_parameter< std::string >::type op(opSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_binary_operation(left_input, right_input, op));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_associative_arithmetic
SEXP apply_delayed_associative_arithmetic(SEXP raw_input, Rcpp::NumericVector val, bool row, std::string op);
RcppExport SEXP _beachmat_apply_delayed_associative_arithmetic(SEXP raw_inputSEXP, SEXP valSEXP, SEXP rowSEXP, SEXP opSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type val(valSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< std::string >::type op(opSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_associative_arithmetic(raw_input, val, row, op));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_nonassociative_arithmetic
SEXP apply_delayed_nonassociative_arithmetic(SEXP raw_input, Rcpp::NumericVector val, bool right, bool row, std::string op);
RcppExport SEXP _beachmat_apply_delayed_nonassociative_arithmetic(SEXP raw_inputSEXP, SEXP valSEXP, SEXP rightSEXP, SEXP rowSEXP, SEXP opSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type val(valSEXP);
    Rcpp::traits::input_parameter< bool >::type right(rightSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< std::string >::type op(opSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_nonassociative_arithmetic(raw_input, val, right, row, op));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_comparison
SEXP apply_delayed_comparison(SEXP raw_input, Rcpp::NumericVector val, bool row, std::string op);
RcppExport SEXP _beachmat_apply_delayed_comparison(SEXP raw_inputSEXP, SEXP valSEXP, SEXP rowSEXP, SEXP opSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type val(valSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< std::string >::type op(opSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_comparison(raw_input, val, row, op));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_boolean
SEXP apply_delayed_boolean(SEXP raw_input, Rcpp::LogicalVector val, bool row, std::string op);
RcppExport SEXP _beachmat_apply_delayed_boolean(SEXP raw_inputSEXP, SEXP valSEXP, SEXP rowSEXP, SEXP opSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type val(valSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    Rcpp::traits::input_parameter< std::string >::type op(opSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_boolean(raw_input, val, row, op));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_boolean_not
SEXP apply_delayed_boolean_not(SEXP raw_input);
RcppExport SEXP _beachmat_apply_delayed_boolean_not(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_boolean_not(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_log
SEXP apply_delayed_log(SEXP raw_input, double base);
RcppExport SEXP _beachmat_apply_delayed_log(SEXP raw_inputSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< double >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_log(raw_input, base));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_log1p
SEXP apply_delayed_log1p(SEXP raw_input);
RcppExport SEXP _beachmat_apply_delayed_log1p(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_log1p(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_abs
SEXP apply_delayed_abs(SEXP raw_input);
RcppExport SEXP _beachmat_apply_delayed_abs(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_abs(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_sqrt
SEXP apply_delayed_sqrt(SEXP raw_input);
RcppExport SEXP _beachmat_apply_delayed_sqrt(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_sqrt(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_ceiling
SEXP apply_delayed_ceiling(SEXP raw_input);
RcppExport SEXP _beachmat_apply_delayed_ceiling(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_ceiling(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_floor
SEXP apply_delayed_floor(SEXP raw_input);
RcppExport SEXP _beachmat_apply_delayed_floor(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_floor(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_round
SEXP apply_delayed_round(SEXP raw_input);
RcppExport SEXP _beachmat_apply_delayed_round(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_round(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_exp
SEXP apply_delayed_exp(SEXP raw_input);
RcppExport SEXP _beachmat_apply_delayed_exp(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_exp(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_subset
SEXP apply_delayed_subset(SEXP raw_input, Rcpp::IntegerVector subset, bool row);
RcppExport SEXP _beachmat_apply_delayed_subset(SEXP raw_inputSEXP, SEXP subsetSEXP, SEXP rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_subset(raw_input, subset, row));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_transpose
SEXP apply_delayed_transpose(SEXP raw_input);
RcppExport SEXP _beachmat_apply_delayed_transpose(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_transpose(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// apply_delayed_bind
SEXP apply_delayed_bind(Rcpp::List input, bool row);
RcppExport SEXP _beachmat_apply_delayed_bind(SEXP inputSEXP, SEXP rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type input(inputSEXP);
    Rcpp::traits::input_parameter< bool >::type row(rowSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_delayed_bind(input, row));
    return rcpp_result_gen;
END_RCPP
}
// initialize_dense_matrix
SEXP initialize_dense_matrix(Rcpp::RObject raw_x, int nrow, int ncol);
RcppExport SEXP _beachmat_initialize_dense_matrix(SEXP raw_xSEXP, SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type raw_x(raw_xSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_dense_matrix(raw_x, nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// fragment_sparse_rows
Rcpp::List fragment_sparse_rows(Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::IntegerVector limits);
RcppExport SEXP _beachmat_fragment_sparse_rows(SEXP iSEXP, SEXP pSEXP, SEXP limitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type limits(limitsSEXP);
    rcpp_result_gen = Rcpp::wrap(fragment_sparse_rows(i, p, limits));
    return rcpp_result_gen;
END_RCPP
}
// sparse_subset_index
Rcpp::IntegerVector sparse_subset_index(Rcpp::IntegerVector starts, Rcpp::IntegerVector newp);
RcppExport SEXP _beachmat_sparse_subset_index(SEXP startsSEXP, SEXP newpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type starts(startsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type newp(newpSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_subset_index(starts, newp));
    return rcpp_result_gen;
END_RCPP
}
// initialize_sparse_matrix
SEXP initialize_sparse_matrix(Rcpp::RObject raw_x, Rcpp::RObject raw_i, Rcpp::RObject raw_p, int nrow, int ncol, bool byrow);
RcppExport SEXP _beachmat_initialize_sparse_matrix(SEXP raw_xSEXP, SEXP raw_iSEXP, SEXP raw_pSEXP, SEXP nrowSEXP, SEXP ncolSEXP, SEXP byrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type raw_x(raw_xSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type raw_i(raw_iSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type raw_p(raw_pSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_sparse_matrix(raw_x, raw_i, raw_p, nrow, ncol, byrow));
    return rcpp_result_gen;
END_RCPP
}
// initialize_SVT_SparseMatrix
SEXP initialize_SVT_SparseMatrix(int nr, int nc, Rcpp::RObject seed);
RcppExport SEXP _beachmat_initialize_SVT_SparseMatrix(SEXP nrSEXP, SEXP ncSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_SVT_SparseMatrix(nr, nc, seed));
    return rcpp_result_gen;
END_RCPP
}
// tatami_dim
Rcpp::IntegerVector tatami_dim(SEXP raw_input);
RcppExport SEXP _beachmat_tatami_dim(SEXP raw_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(tatami_dim(raw_input));
    return rcpp_result_gen;
END_RCPP
}
// tatami_column
Rcpp::NumericVector tatami_column(SEXP raw_input, int i);
RcppExport SEXP _beachmat_tatami_column(SEXP raw_inputSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(tatami_column(raw_input, i));
    return rcpp_result_gen;
END_RCPP
}
// tatami_row
Rcpp::NumericVector tatami_row(SEXP raw_input, int i);
RcppExport SEXP _beachmat_tatami_row(SEXP raw_inputSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(tatami_row(raw_input, i));
    return rcpp_result_gen;
END_RCPP
}
// tatami_row_sums
Rcpp::NumericVector tatami_row_sums(SEXP raw_input, int threads);
RcppExport SEXP _beachmat_tatami_row_sums(SEXP raw_inputSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(tatami_row_sums(raw_input, threads));
    return rcpp_result_gen;
END_RCPP
}
// tatami_column_sums
Rcpp::NumericVector tatami_column_sums(SEXP raw_input, int threads);
RcppExport SEXP _beachmat_tatami_column_sums(SEXP raw_inputSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type raw_input(raw_inputSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(tatami_column_sums(raw_input, threads));
    return rcpp_result_gen;
END_RCPP
}
// initialize_unknown_matrix
SEXP initialize_unknown_matrix(Rcpp::RObject input);
RcppExport SEXP _beachmat_initialize_unknown_matrix(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_unknown_matrix(input));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_beachmat_apply_delayed_binary_operation", (DL_FUNC) &_beachmat_apply_delayed_binary_operation, 3},
    {"_beachmat_apply_delayed_associative_arithmetic", (DL_FUNC) &_beachmat_apply_delayed_associative_arithmetic, 4},
    {"_beachmat_apply_delayed_nonassociative_arithmetic", (DL_FUNC) &_beachmat_apply_delayed_nonassociative_arithmetic, 5},
    {"_beachmat_apply_delayed_comparison", (DL_FUNC) &_beachmat_apply_delayed_comparison, 4},
    {"_beachmat_apply_delayed_boolean", (DL_FUNC) &_beachmat_apply_delayed_boolean, 4},
    {"_beachmat_apply_delayed_boolean_not", (DL_FUNC) &_beachmat_apply_delayed_boolean_not, 1},
    {"_beachmat_apply_delayed_log", (DL_FUNC) &_beachmat_apply_delayed_log, 2},
    {"_beachmat_apply_delayed_log1p", (DL_FUNC) &_beachmat_apply_delayed_log1p, 1},
    {"_beachmat_apply_delayed_abs", (DL_FUNC) &_beachmat_apply_delayed_abs, 1},
    {"_beachmat_apply_delayed_sqrt", (DL_FUNC) &_beachmat_apply_delayed_sqrt, 1},
    {"_beachmat_apply_delayed_ceiling", (DL_FUNC) &_beachmat_apply_delayed_ceiling, 1},
    {"_beachmat_apply_delayed_floor", (DL_FUNC) &_beachmat_apply_delayed_floor, 1},
    {"_beachmat_apply_delayed_round", (DL_FUNC) &_beachmat_apply_delayed_round, 1},
    {"_beachmat_apply_delayed_exp", (DL_FUNC) &_beachmat_apply_delayed_exp, 1},
    {"_beachmat_apply_delayed_subset", (DL_FUNC) &_beachmat_apply_delayed_subset, 3},
    {"_beachmat_apply_delayed_transpose", (DL_FUNC) &_beachmat_apply_delayed_transpose, 1},
    {"_beachmat_apply_delayed_bind", (DL_FUNC) &_beachmat_apply_delayed_bind, 2},
    {"_beachmat_initialize_dense_matrix", (DL_FUNC) &_beachmat_initialize_dense_matrix, 3},
    {"_beachmat_fragment_sparse_rows", (DL_FUNC) &_beachmat_fragment_sparse_rows, 3},
    {"_beachmat_sparse_subset_index", (DL_FUNC) &_beachmat_sparse_subset_index, 2},
    {"_beachmat_initialize_sparse_matrix", (DL_FUNC) &_beachmat_initialize_sparse_matrix, 6},
    {"_beachmat_initialize_SVT_SparseMatrix", (DL_FUNC) &_beachmat_initialize_SVT_SparseMatrix, 3},
    {"_beachmat_tatami_dim", (DL_FUNC) &_beachmat_tatami_dim, 1},
    {"_beachmat_tatami_column", (DL_FUNC) &_beachmat_tatami_column, 2},
    {"_beachmat_tatami_row", (DL_FUNC) &_beachmat_tatami_row, 2},
    {"_beachmat_tatami_row_sums", (DL_FUNC) &_beachmat_tatami_row_sums, 2},
    {"_beachmat_tatami_column_sums", (DL_FUNC) &_beachmat_tatami_column_sums, 2},
    {"_beachmat_initialize_unknown_matrix", (DL_FUNC) &_beachmat_initialize_unknown_matrix, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_beachmat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
