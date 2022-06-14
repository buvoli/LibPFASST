#!/bin/bash

# ---------------------------
# Cori initial setup script (after cloning repository)
# ---------------------------

set -e

git fetch origin
git checkout --track origin/exp-parareal
module load openmpi
make fftw3
make libnpy
make USE_FFT=TRUE USE_FFTW=TRUE

set +e