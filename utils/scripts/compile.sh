#!/bin/bash

gcmDir='/home/glwagner/software/MITgcm'

# Load modules and move to the build directory
module add engaging/intel/2013.1.046
module add netcdf netcdff

cd ./build

# Run genmake with appropriate options, and then make 
$gcmDir/tools/genmake2 \
    -mods=../code/ \
    -optfile=$gcmDir/tools/build_options/linux_amd64_ifort+impi \
    -mpi \
    -rootdir=$gcmDir/

make depend
make -j12
