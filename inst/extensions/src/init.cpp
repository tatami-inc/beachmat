#include "R_ext/Rdynload.h"
#include "integer_access.h"
#include "character_access.h"

#define REGISTER(x) R_RegisterCCallable("morebeach", #x, reinterpret_cast<DL_FUNC>(x))

extern "C" {

void R_init_morebeach(DllInfo *info) {
    // For integers:
    REGISTER(create_integer);
    REGISTER(destroy_integer);
    REGISTER(clone_integer);

    REGISTER(get_dim_integer);
    REGISTER(load_integer);
    REGISTER(load_row2int_integer);
    REGISTER(load_col2int_integer);
    REGISTER(load_row2dbl_integer);
    REGISTER(load_col2dbl_integer);

    REGISTER(load_rows2int_integer);
    REGISTER(load_cols2int_integer);
    REGISTER(load_rows2dbl_integer);
    REGISTER(load_cols2dbl_integer);

    // For characters:
    REGISTER(create_character);
    REGISTER(destroy_character);
    REGISTER(clone_character);

    REGISTER(get_dim_character);
    REGISTER(load_character);
    REGISTER(load_row_character);
    REGISTER(load_col_character);

    REGISTER(load_rows_character);
    REGISTER(load_cols_character);
}

}
