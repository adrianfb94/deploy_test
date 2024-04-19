#!/bin/bash

home=`pwd`
mkdir -p $home/cdo_installed


cd $home/cdo-1.9.1
./configure --with-netcdf=/usr -prefix=$home/cdo_installed
make && make check && make install


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


# #   download, compile and install --> netCDF
# cd $home
# wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.5.0.tar.gz
# tar -xzvf netcdf-4.5.0.tar.gz
# cd netcdf-4.5.0/
# CPPFLAGS=-I~/cdo_install/include LDFLAGS=-L~/cdo_install/lib ./configure -prefix=~/cdo_install CFLAGS=-fPIC
# make && make check && make install

# #   download, compile and install --> jasper
# cd $home
# wget http://www.ece.uvic.ca/~frodo/jasper/software/jasper-1.900.1.zip
# unzip jasper-1.900.1.zip
# cd jasper-1.900.1
# ./configure -prefix=~/cdo_install CFLAGS=-fPIC
# make && make check && make install

# #   download, compile and install --> grib_api
# cd $home
# wget https://software.ecmwf.int/wiki/download/attachments/3473437/grib_api-1.24.0-Source.tar.gz?api=v2 -O grib_api-1.24.0.tar.gz
# tar -xzvf grib_api-1.24.0.tar.gz
# cd grib_api-1.24.0-Source
# ./configure -prefix=~/cdo_install CFLAGS=-fPIC -with-netcdf=~/cdo_install -with-jasper=~/cdo_install
# make && make check && make install

# #   download, compile and install --> cdo
# cd $home
# wget https://code.mpimet.mpg.de/attachments/download/15653/cdo-1.9.1.tar.gz
# tar -xvzf cdo-1.9.1.tar.gz
# cd cdo-1.9.1
# ./configure -prefix=~/cdo_install CFLAGS=-fPIC -with-netcdf=~/cdo_install -with-jasper=~/cdo_install -with-hdf5=~/cdo_install -with-grib_api=~/cdo_install/grib_api-1.24.0-Source
# make && make check && make install

# #   set PATH
# #echo "PATH=\"/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:~/cdo_install/bin\"" > /etc/environment



# # cdo
# cd $home
# ./configure --with-netcdf=/usr -prefix=$home/cdo_installed
# make
# make check
# make install
