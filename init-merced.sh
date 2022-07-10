#!/bin/bash

# -----------------------------------------------------------------------------
# Cori initial setup script, after cloning
#   git clone https://github.com/buvoli/LibPFASST.git
#   cd LibPFASST
#   git fetch origin
#   git checkout --track origin/exp-parareal
# -----------------------------------------------------------------------------

set -e

# check for manually added zip files
if [ ! -f  "fftw-3.3.8.tar.gz" ] ; then
    echo "Manually add fftw-3.3.8.tar.gz (http://fftw.org/fftw-3.3.8.tar.gz) to root directory"
    exit 1
fi

if [ ! -f  "master.zip" ] ; then
    echo "Manually add master.zip (https://github.com/libpfasst/libnpy/archive/refs/heads/master.zip) to root directory"
    exit 1
fi

# overwrite Makefile.external
# - disable wget (wget is firewalled on MERCED cluster -- zip files must be manually copied to LibPFASST root directory)
# - compiles libnpy with ifort
if [ -f "Makefile.external.merced" ]; then
    mv Makefile.external.merced Makefile.external
    # add full path to mkdir (use ; as separater for sed command)
    sed -i 's;MKDIR_PATH = "mkdir";MKDIR_PATH = "/usr/bin/mkdir";g' src/pf_dtype.f90
fi

# make fftw using default compilers
make fftw3

# make remaining using intel openmpi
module load openmpi-2.0/intel
make libnpy
make USE_FFT=TRUE USE_FFTW=TRUE
mv *.mod include/ # added since ifort does not support -J Flag

set +e