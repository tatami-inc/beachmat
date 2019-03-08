pkgname=$1
shift
tmp="$*"
matrices=$(echo $tmp | sed "s/ /|/g")

# Creating the export header.
cat << EOT > exports.h
#include "Rcpp.h"

extern "C" {

EOT

cat *.cpp | egrep "(${matrices})_.*{" | sed -E "s/ [^ ]+,/,/g" | sed -E "s/ [^ ]+\) *\{/\);/" | sed "s/;/;\n/" >> exports.h

echo "}" >> exports.h

# Creating the export file.
cat << EOT > exports.cpp
#include "Rcpp.h"
#include "exports.h"
#include "R_ext/Rdynload.h"

#define REGISTER(x) R_RegisterCCallable("$pkgname", #x, reinterpret_cast<DL_FUNC>(x))

extern "C" {

void R_init_$pkgname(DllInfo *info) {

EOT

cat *.cpp | egrep "(${matrices})_.*{" | sed -E "s/ ?\(.*$//g" | sed -E "s/^.* //" | sed "s/^\(.*\)$/REGISTER(\1);\n/" >> exports.cpp

cat << EOT >> exports.cpp
}

}
EOT
