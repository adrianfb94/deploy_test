#!/bin/bash

# server route
# /opt/render/project/src/child-processes/cdo/
# /opt/render/project/src/child-processes/cdo/cdo-1.9.1/libcdi/app

home=`pwd`

cd $home
tar -xzf zlib-1.2.12.tar.gz
cd zlib-1.2.12
./configure --prefix=$home/cdo_exe
make
make install

rm -r $home/zlib-1.2.12


# cp $home/hdf5-1.14.3.tar.gz $path_cdo_opt/
# cd $path_cdo_opt
# tar -xzf hdf5-1.14.3.tar.gz
# cd hdf5-1.14.3
# ./configure --with-zlib=$path_cdo_opt --prefix=$path_cdo_opt/cdo_exe  CFLAGS=-I$path_cdo_opt/include --disable-shared
# make
# make install


# cp $home/netcdf-c-4.9.2.tar.gz $path_cdo_opt/
# cd $path_cdo_opt
# tar -xzf netcdf-c-4.9.2.tar.gz
# cd Unidata-netcdf-c-d4145f3
# CPPFLAGS=-I$path_cdo_opt/include LDFLAGS=-L$path_cdo_opt/lib ./configure --prefix=$path_cdo_opt/cdo_exe CFLAGS=-I$path_cdo_opt/include --disable-shared
# make
# make install

# cp $home/cdo-1.9.1.tar.gz $path_cdo_opt/
# cd $path_cdo_opt
# tar -xzf cdo-1.9.1.tar.gz
# cd cdo-1.9.1
# ./configure --enable-netcdf4 --enable-zlib --prefix=$path_cdo_opt/cdo_exe --with-netcdf=$path_cdo_opt --with-hdf5=$path_cdo_opt CFLAGS=-I$path_cdo_opt/include --disable-shared
# make
# make install
