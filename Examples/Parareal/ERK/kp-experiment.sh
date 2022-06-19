#!/bin/bash

# =============================================================================
#  Runs the KP experiment from paper
# =============================================================================

RUNNER='mpirun'
FLAGS_ERK='-n 1'
FLAGS_PARA='-n 1'

make DIM=2 # ensure that code is compiled

set -x # (print commands are they run)

# Serial ERK (Reference)
if [ ! -d "dat/kp/serial-ref-P0001/" ] ; then
    $RUNNER $FLAGS_ERK main.2d.exe params/KP-serial-fine.nml nsteps_rk=524288 rho=0 outdir='"kp/serial-ref-"'
fi

# Serial ERK (Fine)
if [ ! -d "dat/kp/serial-fine-P0001/" ] ; then
    $RUNNER $FLAGS_ERK main.2d.exe params/KP-serial-fine.nml nsteps_rk=262144 rho=0.0245436926061703 outdir='"kp/serial-fine-"'
fi

# Parareal ERK NG=1,...,3 and K=0,...,28
for j in {1..3}
do
    for i in {0..28}
    do
        $RUNNER $FLAGS_PARA main.2d.exe params/KP-parareal.nml niters="${i}" nsteps_rk="${j} 32" rho=0.0245436926061703 outdir="\"kp/parareal-ng-${j}-k-${i}-\""
    done
done