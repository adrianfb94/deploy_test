#!/bin/bash

home=`pwd`/child-processes/cdo-1.9.1
cd $home
./configure --with-netcdf=/usr -prefix=$home
make && make check && make install
