#!/bin/bash

set -e
git switch -c exp-parareal
make fftw3
make libnpy
make USE_FFT=TRUE USE_FFTW=TRUE
set +e