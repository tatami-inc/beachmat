#ifndef TATAMI_R_UTILS_HPP
#define TATAMI_R_UTILS_HPP

#include "Rcpp.h"
#include <string>
#include <utility>
#include <stdexcept>
#include <memory>
#include "tatami/tatami.hpp"

namespace tatami_r { 

template<typename Data_, typename Index_>
struct Parsed {
    std::shared_ptr<tatami::Matrix<Data_, Index_> > matrix;

    Rcpp::List contents;
};

inline std::string make_to_string(const Rcpp::RObject& str) {
    Rcpp::StringVector as_str(str);
    if (as_str.size()!=1) { 
        throw std::runtime_error("input RObject should contain a single string");
    }
    return Rcpp::as<std::string>(as_str[0]);
}

inline Rcpp::RObject get_class_object(const Rcpp::RObject& incoming) {
    if (!incoming.isObject()) {
        throw std::runtime_error("object has no 'class' attribute");
    }
    return incoming.attr("class");
}

inline std::string get_class_name(const Rcpp::RObject& incoming) {
    return make_to_string(get_class_object(incoming));
}

inline std::string extract_class_package(const Rcpp::RObject& classname) {
    if (!classname.hasAttribute("package")) {
        throw std::runtime_error("class name has no 'package' attribute");
    }
    return make_to_string(classname.attr("package"));
}

inline std::pair<std::string, std::string> get_class_package(const Rcpp::RObject& incoming) {
    Rcpp::RObject classname=get_class_object(incoming);
    return std::make_pair(make_to_string(classname), extract_class_package(classname));
}

inline std::pair<int, int> parse_dims(Rcpp::RObject dims) {
    if (dims.sexp_type()!=INTSXP) {
        throw std::runtime_error("matrix dimensions should be an integer vector");
    }

    Rcpp::IntegerVector d(dims);
    if (d.size()!=2) {
        throw std::runtime_error("matrix dimensions should be of length 2");
    }

    if (d[0]<0 || d[1]<0) {
        throw std::runtime_error("dimensions should be non-negative");
    }

    return std::make_pair(d[0], d[1]);
}

}

#endif
