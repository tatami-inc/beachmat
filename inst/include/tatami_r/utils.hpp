#ifndef TATAMI_R_UTILS_HPP
#define TATAMI_R_UTILS_HPP

#include "Rcpp.h"
#include <string>
#include <utility>
#include <stdexcept>
#include <memory>
#include "tatami/tatami.hpp"

namespace tatami_r { 

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

}

#endif
