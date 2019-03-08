#include "Rcpp.h"
#include "exports.h"
#include "R_ext/Rdynload.h"

#define REGISTER(x) R_RegisterCCallable("morebeach", #x, reinterpret_cast<DL_FUNC>(x))

extern "C" {

void R_init_morebeach(DllInfo *info) {

REGISTER(AaronMatrix_character_input_create);

REGISTER(AaronMatrix_character_input_destroy);

REGISTER(AaronMatrix_character_input_clone);

REGISTER(AaronMatrix_character_input_dim);

REGISTER(AaronMatrix_character_input_get);

REGISTER(AaronMatrix_character_input_getRow);

REGISTER(AaronMatrix_character_input_getCol);

REGISTER(AaronMatrix_character_input_getConstCol);

REGISTER(AaronMatrix_character_input_getConstColIndexed);

REGISTER(AaronMatrix_character_input_getRows);

REGISTER(AaronMatrix_character_input_getCols);

REGISTER(AaronMatrix_integer_input_create);

REGISTER(AaronMatrix_integer_input_destroy);

REGISTER(AaronMatrix_integer_input_clone);

REGISTER(AaronMatrix_integer_input_dim);

REGISTER(AaronMatrix_integer_input_get);

REGISTER(AaronMatrix_integer_input_getRow_integer);

REGISTER(AaronMatrix_integer_input_getCol_integer);

REGISTER(AaronMatrix_integer_input_getRow_numeric);

REGISTER(AaronMatrix_integer_input_getCol_numeric);

REGISTER(AaronMatrix_integer_input_getConstCol);

REGISTER(AaronMatrix_integer_input_getConstColIndexed);

REGISTER(AaronMatrix_integer_input_getRows_integer);

REGISTER(AaronMatrix_integer_input_getCols_integer);

REGISTER(AaronMatrix_integer_input_getRows_numeric);

REGISTER(AaronMatrix_integer_input_getCols_numeric);

}

}
