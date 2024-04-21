#!/bin/bash

home=`pwd`
mkdir $home/cdo_v1.9
path_cdo=$home/cdo_v1.9

tar -xzvf cdo-1.9.1.tar.gz
cd $home/cdo-1.9.1
./configure --enable-netcdf4  --enable-zlib --prefix=$path_cdo --with-netcdf=/usr
make -j 1
make install
