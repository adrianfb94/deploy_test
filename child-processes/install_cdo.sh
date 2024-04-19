#!/bin/bash

home=`pwd`
# mkdir -p $home/cdo_installed

cd $home/cdo-1.9.1
./configure --with-netcdf=/usr -prefix=$home/cdo-1.9.1
# make
# make check
# make install
