#include "get_all.h"
#include "set_all.h"

extern "C" {

// Row conversion.

SEXP set_row_numeric_to_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_row_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_numeric_to_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_row_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_integer_to_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_row_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_integer_to_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_row_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_logical_to_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_row_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_logical_to_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_row_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

// Column conversion.

SEXP set_col_numeric_to_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_col_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_numeric_to_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_col_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_integer_to_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_col_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_integer_to_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_col_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_logical_to_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_col_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_logical_to_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));
    set_col_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

}
