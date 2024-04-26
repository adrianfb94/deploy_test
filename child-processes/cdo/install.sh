
#!/bin/bash

home=`pwd`

# split -b 1M cdo_dependecies.tar.gz cdo_dependecies.tar.gz.part

cat cdo.part* >> cdo_dependecies.tar.gz

tar -xzvf $home/cdo_dependecies.tar.gz

# path_dep=$home/cdo_full_disable-shared
# path_exe=$home/cdo_exe
# path_cdo=$home/cdo-1.9.1


# # export CPPFLAGS=-I$path_dep/include 
# # export LDFLAGS=-L$path_dep/lib 
# # export CFLAGS=-I$path_dep/include


# tar -xzf $home/netcdf-c-4.9.2.tar.gz
# cd $home/Unidata-netcdf-c-d4145f3
# ./configure --prefix=$path_exe --enable-netcdf4 --disable-hdf5 # --disable-shared 
# make
# make install

# cd $home
# tar -xzf $home/cdo-1.9.1.tar.gz
# cd $home/cdo-1.9.1
# ./configure --prefix=$home/cdo_exe --enable-netcdf4 --with-netcdf=$home/cdo_exe
# make
# make install

# cd $home

rm $home/cdo_dependecies.tar.gz
# rm -r $home/bin
rm -r $home/include
rm -r $home/lib
rm -r $home/share

# rm -r $home/Unidata-netcdf-c-d4145f3
# rm -r $path_dep
# rm -r $path_exe
# rm -r $path_cdo
