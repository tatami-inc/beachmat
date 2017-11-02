#include "utils.h"

namespace beachmat {

std::string make_to_string(const Rcpp::RObject& str) {
    if (str.sexp_type()!=STRSXP || Rf_length(str.get__())!=1) {
        throw std::runtime_error("input RObject should contain a single string");
    }
    return Rcpp::as<std::vector<std::string> >(str).front();
}

void throw_custom_error(const std::string& left, const std::string& classname, const std::string& right) {
    std::stringstream err;
    err << left << classname << right;
    throw std::runtime_error(err.str().c_str());
}

/* Class checks. */

std::string get_class(const Rcpp::RObject& incoming) {
    if (!incoming.isObject()) {
        throw std::runtime_error("object has no class attribute");
    }
    return make_to_string(incoming.attr("class"));
}

Rcpp::RObject get_safe_slot(const Rcpp::RObject& incoming, const std::string& slotname) {
    if (!incoming.hasSlot(slotname)) { 
        std::stringstream err;
        err << "no '" << slotname << "' slot in the " << get_class(incoming) << " object";
        throw std::runtime_error(err.str().c_str()); 
    }
    return incoming.slot(slotname);
}

std::string check_Matrix_class (const Rcpp::RObject& mat, const std::string& expected) {
    std::string mattype=get_class(mat);
    if (!mat.isS4() || mattype.empty() || mattype.substr(1)!=expected) {
        throw_custom_error("matrix should be a *", expected, " object");
    }
    return mattype;
}


/* Type checks */

std::string translate_type(int sexp_type) {
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
            throw std::runtime_error(err.str().c_str());
    }
    return should_be;
}

int reverse_translate_type (const std::string& curtype) {
    if (curtype=="logical") {
        return LGLSXP;
    } else if (curtype=="character") {
        return STRSXP;
    } else if (curtype=="integer") {
        return INTSXP;
    } else if (curtype=="double") {
        return REALSXP;
    }
    std::stringstream err;
    err << "unsupported type'" << curtype << "'";
    throw std::runtime_error(err.str().c_str());
}

int find_sexp_type (const Rcpp::RObject& incoming) {
    if (incoming.isObject()) {
        const std::string classname=get_class(incoming);
        if (classname=="DelayedMatrix") {
            Rcpp::Environment delayenv("package:DelayedArray");
            Rcpp::Function typefun=delayenv["type"];
            std::string curtype=Rcpp::as<std::string>(typefun(incoming));
            return reverse_translate_type(curtype);
            
        } else if (classname=="HDF5Matrix") {
            Rcpp::RObject h5seed=get_safe_slot(incoming, "seed");
            Rcpp::RObject first_val=get_safe_slot(h5seed, "first_val");
            return first_val.sexp_type();

        } else if (classname=="RleMatrix") {
            Rcpp::RObject rleseed=get_safe_slot(incoming, "seed");
            std::string rclass=get_class(rleseed);
            if (rclass=="SolidRleArraySeed") {
                Rcpp::RObject rle=get_safe_slot(rleseed, "rle");
                return get_safe_slot(rle, "values").sexp_type();
            } else if (rclass=="ChunkedRleArraySeed") {
                std::string curtype=Rcpp::as<std::string>(get_safe_slot(rleseed, "type"));
                return reverse_translate_type(curtype);               
            }

        } else if (classname.length()==9 && classname.substr(3)=="Matrix") {
            if (classname[0]=='d') {
                return REALSXP;
            } else if (classname[0]=='l') {
                return LGLSXP;
            }

        }
        throw_custom_error("unknown SEXP type for ", classname, " object");
    }
    return incoming.sexp_type();
}

/* DelayedArray utilities. */

bool is_pristine_delayed_array(const Rcpp::RObject& in) {
    const Rcpp::Environment env=Rcpp::Environment::namespace_env("DelayedArray");
    Rcpp::Function fun=env["is_pristine"];
    Rcpp::LogicalVector out=fun(in);
    if (out.size()!=1) {
        throw std::runtime_error("pristine check should return a logical scalar");
    }
    return out[0];
}

Rcpp::RObject realize_delayed_array(const Rcpp::RObject& in) {
    if (!in.isS4() || get_class(in)!="DelayedMatrix") { 
        throw std::runtime_error("object should be a DelayedMatrix");
    }

    // Checking if it is only Delayed because of dimnames on a HDF5Matrix, for which we just construct the seed directly.
    Rcpp::RObject seed(get_safe_slot(in, "seed"));
    if (get_class(seed)=="HDF5ArraySeed") { 
        bool unchanged=true;
        Rcpp::IntegerVector rawdim(get_safe_slot(seed, "dim"));
        if (rawdim.size()!=2) { 
            throw std::runtime_error("seed dimensions should be a vector of length 2");
        }

        // Indices are consecutive and purely increasing, or NULL. 
        Rcpp::List indices(get_safe_slot(in, "index"));
        if (indices.size()!=2) { 
            throw std::runtime_error("indices should be a list of two integer vectors");
        }

        for (size_t d=0; d<2; ++d) { 
            Rcpp::RObject curobj(indices[d]);
            if (curobj.isNULL()) {
                continue;
            }
            
            Rcpp::IntegerVector curvec(curobj);
            if (curvec.size()!=rawdim[d]) { 
                unchanged=false;
                break;
            }

            int count=1;
            for (auto idex : curvec) { 
                if (idex!=count) {
                    unchanged=false;
                    break;
                }
                ++count;
            }
            if (!unchanged) {
                break;
            }
        }
       
        // Delayed operations are empty. 
        Rcpp::List delops=in.slot("delayed_ops");
        if (delops.size()) { 
            unchanged=false;
        }

        // No transposition is performed.
        Rcpp::LogicalVector trans=in.slot("is_transposed");
        if (trans.size()!=1) { 
            throw std::runtime_error("transposition specification should be a logical scalar");
        }
        if (trans[0]) { 
            unchanged=false;
        }
        
        if (unchanged) {
            std::string matclass="HDF5Matrix";
            Rcpp::S4 h5mat(matclass);
            if (!h5mat.hasSlot("seed")) {
                throw_custom_error("missing 'seed' slot in ", matclass, " object");
            }
            h5mat.slot("seed") = in.slot("seed");
            return Rcpp::RObject(h5mat);
        }
    }

    Rcpp::Environment delayenv("package:DelayedArray");
    Rcpp::Function realfun=delayenv["realize"];
    return realfun(in);
}

}
