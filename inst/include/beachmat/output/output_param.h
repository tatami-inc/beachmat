#ifndef BEACHMAT_OUTPUT_PARAM_H
#define BEACHMAT_OUTPUT_PARAM_H

#include "Rcpp.h"

#include "../utils/utils.h"

namespace beachmat {

class output_param {
public:
    output_param() {} // TODO: add option to fetch getRealizationBackend().

    output_param(const std::string& m, const std::string& p) : cls(m), pkg(p) {}

    output_param(Rcpp::RObject in) {
        if (!in.isS4()) {
            return;
        }

        auto classinfo=get_class_package(in);
        cls=classinfo.first;
        pkg=classinfo.second;
        
        return;
    }

    bool is_external_available(const std::string& type) const {
        // Skipping packages that we know won't support this.
        if (pkg=="" || pkg=="Matrix" || pkg=="base" || pkg=="DelayedArray") { return false; }

        // Incidentally ensures that 'pkg' is loaded so its registered symbols
        // are on the search path - otherwise this just crashes.
        Rcpp::Environment pkgenv=Rcpp::Environment::namespace_env(pkg);

        std::stringstream symbolic;
        symbolic << cls << "_" << type << "_output";
        Rcpp::RObject out=pkgenv.get(symbolic.str());

        if (out.isNULL()) {
            return false;
        }
        
        Rcpp::LogicalVector flag(out);
        if (flag.size()!=1) {
            std::stringstream msg;
            msg << "invalid output specifier for " << symbolic.str();
            throw std::runtime_error(msg.str());
        }

        return flag[0];
    }

    std::string get_class() const { return cls; }
    std::string get_package() const { return pkg; }
private:
    std::string cls="matrix";
    std::string pkg="base";
};

}

#endif
