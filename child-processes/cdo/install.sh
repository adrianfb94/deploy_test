#!/bin/bash

home=`pwd`

# split -b 1M cdo_dependecies.tar.gz cdo_dependecies.tar.gz.part

cat cdo_dependecies.tar.gz.part* >> cdo_dependecies.tar.gz

tar -xzvf $home/cdo_dependecies.tar.gz

path_dep=$home/cdo_full_disable-shared
path_exe=$home/cdo_exe
path_cdo=$home/cdo-1.9.1


echo $path_dep
echo $path_exe
echo $path_cdo

rm $home/cdo_dependecies.tar.gz
rm -r $path_dep
rm -r $path_exe
rm -r $path_cdo

# cd $path_cdo
# ./configure --enable-netcdf4 --enable-zlib --prefix=$path_exe --with-netcdf=$path_dep --with-hdf5=$path_dep CFLAGS=-I$path_dep/include
# make
# make install

