#ifndef EXTERNAL_READER_H
#define EXTERNAL_READER_H

#include "Rcpp.h"

#include "dim_checker.h"

#include <string>
#include <sstream>

namespace beachmat {

// Assistant function to define names.
inline std::string get_external_name(const char* matclass, const char* type, const char* mode, const char* fname, const char* intype=NULL) {
    std::stringstream fname;
    fname << matclass << "_" << type << "_" << mode << "_" << fname;
    if (intype) {
        fname << "_" << intype;
    }
    return fname.str();
}
 
// Carefully copied external pointer.
class external_ptr {
private:
    void* ptr=NULL;
    void * (*clone) (void *)=NULL;
    void (*destroy) (void *)=NULL;
    void * safe_clone() const { return (ptr!=NULL ? clone(ptr) : NULL); }
    void self_destruct() const { 
        if (ptr!=NULL) { destroy(ptr); } 
        return;
    }
public:
    external_ptr() = default;
    ~external_ptr() {
        self_destruct();
        return;
    }

    external_ptr(SEXP in, const char* pkg, const char* matclass, const char* type) { // input constructor
        auto clone_name=get_external_name(matclass, type, "input", "clone");
        clone=reinterpret_cast<void * (*)(void *)>(R_GetCCallable(pkg, clone_name.c_str()));

        auto destroy_name=get_external_name(matclass, type, "input", "destroy");
        destroy=reinterpret_cast<void (*)(void *)>(R_GetCCallable(pkg, destroy_name.c_str()));

        auto create_name=get_external_name(matclass, type, "input", "create");
        auto create=reinterpret_cast<void * (*)(SEXP)>(R_GetCCallable(pkg, create_name.c_str()));
        ptr=create(in);
        return;
    }

    external_ptr(size_t nr, size_t nc, const char* pkg, const char* matclass, const char* type) { // output constructor
        auto clone_name=get_external_name(matclass, type, "output", "clone");
        clone=reinterpret_cast<void * (*)(void *)>(R_GetCCallable(pkg, clone_name.c_str()));

        auto destroy_name=get_external_name(matclass, type, "output", "destroy");
        destroy=reinterpret_cast<void (*)(void *)>(R_GetCCallable(pkg, destroy_name.c_str()));

        auto create_name=get_external_name(matclass, type, "output", "create");
        auto create=reinterpret_cast<void * (*)()>(R_GetCCallable(pkg, create_name.c_str()));
        ptr=create(nr, nc);
        return;
    }

    external_ptr(const external_ptr& other) : ptr(other.safe_clone()), clone(other.clone), destroy(other.destroy) {}
    external_ptr& operator=(const external_ptr& other) {
        self_destruct(); // avoid memory leak.
        ptr=other.safe_clone();
        clone=other.clone;
        destroy=other.destroy;
        return *this;
    }
    external_ptr(external_ptr&& other) : ptr(other.ptr), clone(other.clone), destroy(other.destroy) {
        other.ptr=NULL; // avoid double destruction.
        other.clone=NULL;
        other.destroy=NULL;
        return;
    }
    external_ptr& operator=(external_ptr&& other) {
        self_destruct(); // avoid memory leak.
        ptr=other.ptr;
        clone=other.clone;
        destroy=other.destroy;

        other.ptr=NULL; // avoid double destruction.
        other.clone=NULL;
        other.destroy=NULL;
        return *this;
    }
    void* get() const { return ptr; }
};

}

#endif
