cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --prefix=$PREFIX
