#!/bin/bash

# -----------------------------------------------------------------------------
# Local initial setup script, after cloning
#   git clone https://github.com/buvoli/LibPFASST.git
#   git fetch origin
#   git checkout --track origin/exp-parareal
# -----------------------------------------------------------------------------

set -e
make fftw3
make libnpy
make USE_FFT=TRUE USE_FFTW=TRUE
set +e