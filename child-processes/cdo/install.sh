#!/bin/bash

home=`pwd`



# split -b 1M cdo_dependecies.tar.gz cdo_dependecies.tar.gz.part

cat cdo_dependecies.tar.gz.part* >> cdo_dependecies.tar.gz

tar -xzf $home/cdo_dependecies.tar.gz

path_dep=$home/cdo_full_disable-shared
# path_dep=/opt/render/project/src/.venv/lib/python3.11/site-packages/netCDF4.libs
path_exe=$home/cdo_exe
path_cdo=$home/cdo-1.9.1


cd $path_cdo
# CPPFLAGS=-I$path_dep/include LDFLAGS=-L$path_dep/lib CFLAGS=-I$path_dep/include ./configure --prefix=$path_exe --with-netcdf=$path_dep # --with-hdf5=$path_dep
./configure --prefix=$path_exe --with-netcdf=$path_dep # --with-hdf5=$path_dep
# ./configure --enable-netcdf4 --enable-zlib --prefix=$path_exe --with-netcdf=$path_dep --with-hdf5=$path_dep

# make
# make install

rm $home/cdo_dependecies.tar.gz
rm -r $path_dep
rm -r $path_exe
rm -r $path_cdo
