#ifndef BEACHMAT_HDF5_READER_H
#define BEACHMAT_HDF5_READER_H

#include "beachmat.h"
#include "dim_checker.h"
#include "HDF5_utils.h"

namespace beachmat {

/*** Class definition ***/

template<typename T, int RTYPE>
class HDF5_reader : public dim_checker {
public:
    HDF5_reader(const Rcpp::RObject&);
    ~HDF5_reader();

    template<typename X>
    void extract_row(size_t, X*, const H5::DataType&, size_t, size_t);

    template<typename X>
    void extract_col(size_t, X*, const H5::DataType&, size_t, size_t);
    
    template<typename X>
    void extract_one(size_t, size_t, X*, const H5::DataType&);  

    template<typename X>
    void extract_rows(Rcpp::IntegerVector::iterator, size_t, X*, const H5::DataType&, size_t, size_t);

    template<typename X>
    void extract_cols(Rcpp::IntegerVector::iterator, size_t, X*, const H5::DataType&, size_t, size_t);

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type() const;

    H5::DataType get_datatype() const;
protected:
    Rcpp::RObject original;
    std::string filename, dataname;

    H5::H5File hfile;
    H5::DataSet hdata;
    HDF5_selector hselect;

    bool onrow, oncol;
    bool rowokay, colokay;
    bool largerrow, largercol;
    H5::FileAccPropList rowlist, collist;
};

/*** Constructor definition ***/

template<typename T, int RTYPE>
HDF5_reader<T, RTYPE>::HDF5_reader(const Rcpp::RObject& incoming) : original(incoming),
        rowlist(H5::FileAccPropList::DEFAULT.getId()), collist(H5::FileAccPropList::DEFAULT.getId()) {

    std::string ctype=get_class(incoming);
    if (!incoming.isS4() || ctype!="HDF5Matrix") {
        throw std::runtime_error("matrix should be a HDF5Matrix or DelayedMatrix object");
    }

    const Rcpp::RObject& h5_seed=get_safe_slot(incoming, "seed");
    std::string stype=get_class(h5_seed);
    if (!h5_seed.isS4() || stype!="HDF5ArraySeed") {
        throw_custom_error("'seed' slot in a ", ctype, " object should be a HDF5ArraySeed object");
    }

    // Checking first value.
    const Rcpp::RObject& firstval=get_safe_slot(h5_seed, "first_val");
    if (firstval.sexp_type()!=RTYPE) { 
        std::stringstream err;
        err << "'first_val' slot in a " << get_class(h5_seed) << " object should be " << translate_type(RTYPE);
        throw std::runtime_error(err.str().c_str());
    }

    // Checking dimensions.
    if (!h5_seed.hasAttribute("dim")) { 
        throw_custom_error("", stype, " object should have 'dim' attribute"); 
    }
    this->fill_dims(h5_seed.attr("dim"));
    const size_t& NC=this->ncol;
    const size_t& NR=this->nrow;

    // Checking names.
    try {
        filename=make_to_string(get_safe_slot(h5_seed, "filepath"));
    } catch (...) { 
        throw_custom_error("'filepath' slot in a ", stype, " object should be a string");
    }
    try {
        dataname=make_to_string(get_safe_slot(h5_seed, "name"));
    } catch (...) { 
        throw_custom_error("'name' slot in a ", stype, " object should be a string");
    }
    
    // Setting up the HDF5 accessors.
    hfile.openFile(filename.c_str(), H5F_ACC_RDONLY);
    hdata = hfile.openDataSet(dataname.c_str());

    H5::DataSpace hspace = hdata.getSpace();
    if (hspace.getSimpleExtentNdims()!=2) {
        throw std::runtime_error("data in HDF5 file is not a two-dimensional array");
    }

    hsize_t dims_out[2];
    hspace.getSimpleExtentDims(dims_out, NULL);
    if (dims_out[1]!=NR || dims_out[0]!=NC) { 
        throw_custom_error("dimensions in HDF5 file do not equal dimensions in the ", ctype, " object");
    }

    hselect.set_dims(NR, NC);

    // Checking the type.
    auto curtype=hdata.getTypeClass();
    switch (RTYPE) {
        case REALSXP:
            if (curtype!=H5T_FLOAT) { 
                throw std::runtime_error("data type in HDF5 file is not double");
            }
            break;
        case INTSXP: 
            if (curtype!=H5T_INTEGER) { 
                throw std::runtime_error("data type in HDF5 file is not integer");
            }
            break;
        case LGLSXP:
            if (curtype!=H5T_INTEGER) { 
                throw std::runtime_error("data type in HDF5 file is not logical");
            }
            break;
        case STRSXP:
            if (curtype!=H5T_STRING) { 
                throw std::runtime_error("data type in HDF5 file is not character");
            }
            break;
    }

    // Setting the chunk cache parameters.
    calc_HDF5_chunk_cache_settings(this->nrow, this->ncol, hdata.getCreatePlist(), this->get_datatype(),
            onrow, oncol, rowokay, colokay, largerrow, largercol, rowlist, collist);
    return;
}

template<typename T, int RTYPE>
HDF5_reader<T, RTYPE>::~HDF5_reader() {}

/*** Basic getter methods ***/

template<typename T, int RTYPE>
template<typename X>
void HDF5_reader<T, RTYPE>::extract_row(size_t r, X* out, const H5::DataType& HDT, size_t first, size_t last) { 
    check_rowargs(r, first, last);
    reopen_HDF5_file_by_dim(filename, dataname, 
            hfile, hdata, H5F_ACC_RDONLY, rowlist, 
            onrow, oncol, largercol, rowokay);
    hselect.select_row(r, first, last);
    hdata.read(out, HDT, hselect.get_row_space(), hselect.get_mat_space());
    return;
}

template<typename T, int RTYPE>
template<typename X>
void HDF5_reader<T, RTYPE>::extract_col(size_t c, X* out, const H5::DataType& HDT, size_t first, size_t last) { 
    check_colargs(c, first, last);
    reopen_HDF5_file_by_dim(filename, dataname, 
            hfile, hdata, H5F_ACC_RDONLY, collist, 
            oncol, onrow, largerrow, colokay);
    hselect.select_col(c, first, last);
    hdata.read(out, HDT, hselect.get_col_space(), hselect.get_mat_space());
    return;
}
    
template<typename T, int RTYPE>
template<typename X>
void HDF5_reader<T, RTYPE>::extract_one(size_t r, size_t c, X* out, const H5::DataType& HDT) { 
    check_oneargs(r, c);
    hselect.select_one(r, c);
    hdata.read(out, HDT, hselect.get_one_space(), hselect.get_mat_space());
    return;
}

/*** Multi getter methods ***/

template<typename T, int RTYPE>
template<typename X>
void HDF5_reader<T, RTYPE>::extract_rows(Rcpp::IntegerVector::iterator cIt, size_t n, X* out, const H5::DataType& HDT, size_t first, size_t last) { 
    check_rowargs(0, first, last);
    check_row_indices(cIt, n);

    reopen_HDF5_file_by_dim(filename, dataname, 
            hfile, hdata, H5F_ACC_RDONLY, rowlist, 
            onrow, oncol, largercol, rowokay);

    hselect.clear_mat_space();
    for (size_t i=0; i<n; ++i, ++cIt) {
        hselect.select_row(*cIt, first, last, H5S_SELECT_OR);
    }

    hsize_t custom_dim[2];
    custom_dim[0]=last-first;
    custom_dim[1]=n;
    H5::DataSpace custom_space(2, custom_dim);
    hdata.read(out, HDT, custom_space, hselect.get_mat_space()); 
    return;
}

template<typename T, int RTYPE>
template<typename X>
void HDF5_reader<T, RTYPE>::extract_cols(Rcpp::IntegerVector::iterator cIt, size_t n, X* out, const H5::DataType& HDT, size_t first, size_t last) {
    check_colargs(0, first, last);
    check_col_indices(cIt, n);

    reopen_HDF5_file_by_dim(filename, dataname, 
            hfile, hdata, H5F_ACC_RDONLY, collist, 
            oncol, onrow, largerrow, colokay);

    hselect.clear_mat_space();
    for (size_t i=0; i<n; ++i, ++cIt) {
        hselect.select_col(*cIt, first, last, H5S_SELECT_OR);
    }

    hsize_t custom_dim[2];
    custom_dim[0]=n;
    custom_dim[1]=last-first;
    H5::DataSpace custom_space(2, custom_dim);
    hdata.read(out, HDT, custom_space, hselect.get_mat_space()); 
    return;
}

/*** Other methods ***/

template<typename T, int RTYPE>
H5::DataType HDF5_reader<T, RTYPE>::get_datatype() const {
    return hdata.getDataType();
}

template<typename T, int RTYPE>
Rcpp::RObject HDF5_reader<T, RTYPE>::yield() const {
    return original;
}

template<typename T, int RTYPE>
matrix_type HDF5_reader<T, RTYPE>::get_matrix_type() const {
    return HDF5;
}

}

#endif
