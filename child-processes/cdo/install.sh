#!/bin/bash

# server route
# /opt/render/project/src/child-processes/cdo/
# /opt/render/project/src/child-processes/cdo/cdo-1.9.1/libcdi/app

home=`pwd`

# cd $home
# tar -xzf zlib-1.2.12.tar.gz
# cd zlib-1.2.12
# ./configure --prefix=$home/cdo_exe
# make
# make install

# rm -r $home/zlib-1.2.12

cd $home
tar -xzf hdf5-1.14.3.tar.gz
cd hdf5-1.14.3
./configure --with-zlib=$home/cdo_exe --prefix=$home/cdo_exe 
make
make install

rm -r $home/hdf5-1.14.3


cd $home
tar -xzf netcdf-c-4.9.2.tar.gz
cd Unidata-netcdf-c-d4145f3
./configure --prefix=$home/cdo_exe
make
make install

rm -r $home/Unidata-netcdf-c-d4145f3



cd $home
tar -xzf cdo-1.9.1.tar.gz
cd cdo-1.9.1
./configure --enable-netcdf4 --enable-zlib --prefix=$home/cdo_exe --with-netcdf=$home/cdo_exe --with-hdf5=$home/cdo_exe
make
make install
