#!/bin/bash

# -----------------------------------------------------------------------------
# Cori initial setup script, after cloning
#   git clone https://github.com/buvoli/LibPFASST.git
#   cd LibPFASST
#   git fetch origin
#   git checkout --track origin/exp-parareal
# -----------------------------------------------------------------------------

set -e

module load openmpi

mv Makefile.defaults.cori Makefile.defaults # disable sweepers
mv Makefile.external.cori Makefile.external # compiles libnpy with ifort

make fftw3
make libnpy
make USE_FFT=TRUE USE_FFTW=TRUE
mv *.mod include/ # added since ifort does not support -J Flag

set +e