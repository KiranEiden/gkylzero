#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/OpenBLAS-0.3.23
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading OpenBLAS .."
    # delete old checkout and builds
    rm -rf OpenBLAS-*
    curl -L https://github.com/xianyi/OpenBLAS/releases/download/v0.3.23/OpenBLAS-0.3.23.tar.gz > OpenBLAS-0.3.23.tar.gz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building OpenBLAS .."
    gunzip -f OpenBLAS-0.3.23.tar.gz
    tar xvf OpenBLAS-0.3.23.tar
    cd OpenBLAS-0.3.23
    make USE_OPENMP=0 NUM_THREADS=1 FC=$FC -j 32
    make USE_OPENMP=0 NUM_THREADS=1 install PREFIX=$PREFIX -j 32

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/OpenBLAS
fi


