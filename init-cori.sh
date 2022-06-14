#!/bin/bash

# -----------------------------------------------------------------------------
# Cori initial setup script, after cloning
#   git clone https://github.com/buvoli/LibPFASST.git
#   git fetch origin
#   git checkout --track origin/exp-parareal
# -----------------------------------------------------------------------------

set -e

module load openmpi
make fftw3
make libnpy
make USE_FFT=TRUE USE_FFTW=TRUE

set +e