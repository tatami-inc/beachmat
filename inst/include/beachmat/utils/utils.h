#ifndef BEACHMAT_UTILS_H
#define BEACHMAT_UTILS_H

#include "Rcpp.h"

#include <string>
#include <sstream>
#include <utility>
#include <tuple>
#include <stdexcept>

namespace beachmat { 

// Typedef for the indexing tuple.

template<class V>
using const_col_indexed_info=std::tuple<size_t, Rcpp::IntegerVector::iterator, typename V::iterator>;

// Definition of a copyable vector, to use in classes.

template <class V> 
struct copyable_holder {
    copyable_holder(size_t n=0) : vec(n) {}
    ~copyable_holder() = default;
    copyable_holder(const copyable_holder& other) : vec(Rcpp::clone(other.vec)) {}
    copyable_holder& operator=(const copyable_holder& other) { 
        vec=Rcpp::clone(other.vec); 
        return *this;
    }
    copyable_holder(copyable_holder&&) = default;
    copyable_holder& operator=(copyable_holder&&) = default;
    V vec;
};

/* String-related helper functions */

inline std::string make_to_string(const Rcpp::RObject& str) {
    Rcpp::StringVector as_str(str);
    if (as_str.size()!=1) { 
        throw std::runtime_error("input RObject should contain a single string");
    }
    return Rcpp::as<std::string>(as_str[0]);
}

template <class L, class R>
std::string combine_strings(const L& left, const R& right) {
    std::stringstream err;
    err << left << right;
    return err.str();    
}

inline void throw_custom_error(const std::string& left, const std::string& classname, const std::string& right) {
    std::stringstream err;
    err << left << classname << right;
    throw std::runtime_error(err.str());
}

/* Class checks. */

inline Rcpp::RObject get_class_object(const Rcpp::RObject& incoming) {
    if (!incoming.isObject()) {
        throw std::runtime_error("object has no 'class' attribute");
    }
    return incoming.attr("class");
}

inline std::string get_class(const Rcpp::RObject& incoming) {
    return make_to_string(get_class_object(incoming));
}

inline std::pair<std::string, std::string> get_class_package(const Rcpp::RObject& incoming) {
    Rcpp::RObject classname=get_class_object(incoming);
    if (!classname.hasAttribute("package")) {
        throw std::runtime_error("class name has no 'package' attribute");
    }
    return std::make_pair(make_to_string(classname), make_to_string(classname.attr("package")));
}

inline Rcpp::RObject get_safe_slot(const Rcpp::RObject& incoming, const std::string& slotname) {
    if (!incoming.hasSlot(slotname)) { 
        std::stringstream err;
        err << "no '" << slotname << "' slot in the " << get_class(incoming) << " object";
        throw std::runtime_error(err.str()); 
    }
    return incoming.slot(slotname);
}

/* Type checks */

inline std::string translate_type(int sexp_type) {
    std::string should_be;
    switch(sexp_type) {
        case REALSXP:
            should_be="double";
            break;
        case INTSXP:
            should_be="integer";
            break;
        case LGLSXP:
            should_be="logical";
            break;
        case STRSXP:
            should_be="character";
            break;
        default:
            std::stringstream err;
            err << "unsupported sexptype '" << sexp_type << "'";
            throw std::runtime_error(err.str());
    }
    return should_be;
}

inline int find_sexp_type (const Rcpp::RObject& incoming) {
    if (!incoming.isObject()) {
        return incoming.sexp_type();
    }

    const auto classinfo=get_class_package(incoming);
    const std::string& classname=classinfo.first;
    const std::string& classpkg=classinfo.second;
    
    if (classpkg=="Matrix" && classname.length()==9 && classname.substr(3)=="Matrix") {
        if (classname[0]=='d') {
            return REALSXP;
        } else if (classname[0]=='l') {
            return LGLSXP;
        }

    } else {
        Rcpp::Environment delayenv=Rcpp::Environment::namespace_env("DelayedArray");
        Rcpp::Function typefun=delayenv["type"];
        std::string curtype=Rcpp::as<std::string>(typefun(incoming));
        if (curtype=="logical") {
            return LGLSXP;
        } else if (curtype=="character") {
            return STRSXP;
        } else if (curtype=="integer") {
            return INTSXP;
        } else if (curtype=="double") {
            return REALSXP;
        }
    } 
    throw_custom_error("unknown SEXP type for ", classname, " object");
}

}

#endif
