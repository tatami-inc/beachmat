#include "set_errors.h"
#include "get_errors.h"

extern "C" { 

SEXP set_errors_integer (SEXP in, SEXP mode, SEXP reget) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));

    if (Rcpp::LogicalVector(reget)[0]) {
        get_errors<Rcpp::IntegerVector>(optr.get(), mode);
    } else {
        set_errors<Rcpp::IntegerVector>(optr.get(), mode);
    }
    return optr->yield();
    END_RCPP
}

SEXP set_errors_logical (SEXP in, SEXP mode, SEXP reget) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));

    if (Rcpp::LogicalVector(reget)[0]) {
        get_errors<Rcpp::LogicalVector>(optr.get(), mode);
    } else {
        set_errors<Rcpp::LogicalVector>(optr.get(), mode);
    }
    return optr->yield();
    END_RCPP
}

SEXP set_errors_numeric (SEXP in, SEXP mode, SEXP reget) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));

    if (Rcpp::LogicalVector(reget)[0]) {
        get_errors<Rcpp::NumericVector>(optr.get(), mode);
    } else {
        set_errors<Rcpp::NumericVector>(optr.get(), mode);
    }
    return optr->yield();
    END_RCPP
}

SEXP set_errors_character (SEXP in, SEXP mode, SEXP reget) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), beachmat::output_param(in, false, true));

    if (Rcpp::LogicalVector(reget)[0]) {
        get_errors<Rcpp::StringVector>(optr.get(), mode);
    } else {
        set_errors<Rcpp::StringVector>(optr.get(), mode);
    }
    return optr->yield();
    END_RCPP
}

}
