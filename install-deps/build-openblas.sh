#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/OpenBLAS-0.3.15
# Location where dependency sources will be downloaded
DEP_SOURCES=$HOME/gkylsoft/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

# delete old checkout and builds
rm -rf cln-*

curl -L https://github.com/xianyi/OpenBLAS/releases/download/v0.3.15/OpenBLAS-0.3.15.tar.gz > OpenBLAS-0.3.15.tar.gz
gunzip -f OpenBLAS-0.3.15.tar.gz
tar xvf OpenBLAS-0.3.15.tar
cd OpenBLAS-0.3.15
make USE_OPENMP=0 FC=$FC -j
make install PREFIX=$PREFIX -j

# soft-link 
ln -sfn $PREFIX $GKYLSOFT/OpenBLAS


