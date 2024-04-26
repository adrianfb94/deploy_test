#!/bin/bash

home=`pwd`

# split -b 1M cdo_dependecies.tar.gz cdo_dependecies.tar.gz.part

cat cdo_dependecies.tar.gz.part* >> cdo_dependecies.tar.gz

tar -xzf $home/cdo_dependecies.tar.gz

path_dep=$home/cdo_full_disable-shared
path_exe=$home/cdo_exe
path_cdo=$home/cdo-1.9.1


# # export CPPFLAGS=-I$path_dep/include 
# # export LDFLAGS=-L$path_dep/lib 
# # export CFLAGS=-I$path_dep/include


tar -xzf $home/netcdf-c-4.9.2.tar.gz
cd $home/Unidata-netcdf-c-d4145f3
./configure --prefix=$path_exe --enable-netcdf4 --disable-hdf5 # --disable-shared 
make
make install

# cd $path_cdo
# # CPPFLAGS=-I$path_dep/include LDFLAGS=-L$path_dep/lib CFLAGS=-I$path_dep/include ./configure --prefix=$path_exe --with-netcdf=$path_dep # --with-hdf5=$path_dep
# ./configure --prefix=$path_exe --enable-netcdf4 
# # ./configure --enable-netcdf4 --enable-zlib --prefix=$path_exe --with-netcdf=$path_dep --with-hdf5=$path_dep

# # make
# # make install



rm $home/cdo_dependecies.tar.gz
rm -r $home/Unidata-netcdf-c-d4145f3
rm -r $path_dep
# rm -r $path_exe
rm -r $path_cdo

