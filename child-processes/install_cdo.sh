#!/bin/bash

home=`pwd`

# cd $home/zlib-1.2.8
# ./configure -prefix=$home/cdo_installed
# make
# make check
# make install


# cd $home/hdf5-1.8.13
# ./configure -with-zlib=$home/cdo_installed -prefix=$home/cdo_installed CFLAGS=-fPIC
# make
# make check
# make install


# cd $home/netcdf-c-4.5.0
# CPPFLAGS=-I$home/cdo_installed/include LDFLAGS=-L$home/cdo_installed/lib ./configure -prefix=$home/cdo_installed CFLAGS=-fPIC
# make && make check && make install

# cd $home/jasper-1.900.1
# ./configure -prefix=$home/cdo_installed CFLAGS=-fPIC
# make && make check && make install

# cd $home/grib_api-1.24.0-Source
# ./configure -prefix=$home/cdo_installed CFLAGS=-fPIC -with-netcdf=$home/cdo_installed -with-jasper=$home/cdo_installed
# make && make check && make install


# cd $home/cdo-1.9.1
# ./configure -prefix=$home/cdo_installed CFLAGS=-fPIC -with-netcdf=$home/cdo_installed -with-jasper=$home/cdo_installed -with-hdf5=$home/cdo_installed -with-grib_api=$home/grib_api-1.24.0-Source
# make && make check && make install

# #   set PATH
# #echo "PATH=\"/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:~/cdo_install/bin\"" > /etc/environment


cd $home/cdo-1.9.1
./configure -prefix=$home/cdo_installed CFLAGS=-fPIC -with-netcdf=/usr -with-jasper=/usr -with-hdf5=/usr
make && make check && make install
