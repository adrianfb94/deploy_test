#!/bin/bash

tar -xzf cdo-1.9.9.tar.gz
cd cdo-1.9.9
./configure 
make 
make install
