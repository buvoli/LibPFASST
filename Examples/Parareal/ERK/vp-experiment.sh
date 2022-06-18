#!/bin/bash

# =============================================================================
#  Runs the VP experiment from paper
# =============================================================================

RUNNER='mpirun'
FLAGS_ERK='-n 1'
FLAGS_PARA='-n 1'

make DIM=2 # ensure that code is compiled

set -x # (print commands are they run)

# Serial ERK (Reference)
if [ ! -d "dat/vp/serial-ref-P0001/" ] ; then
    $RUNNER $FLAGS_ERK main.2d.exe params/VP-serial-fine.nml nsteps_rk=262144 rho=0 outdir='"vp/serial-ref-"'
fi

# Serial ERK (Fine)
if [ ! -d "dat/vp/serial-fine-P0001/" ] ; then
    $RUNNER $FLAGS_ERK main.2d.exe params/VP-serial-fine.nml nsteps_rk=131072 rho=0.0245436926061703 outdir='"vp/serial-fine-"'
fi

# Parareal ERK NG=1,...,3 and K=0,...,12
for j in {1..3}
do
    for i in {0..12}
    do
        $RUNNER $FLAGS_PARA main.2d.exe params/VP-parareal.nml niters="${i}" nsteps_rk="${j} 64" rho=0.0245436926061703 outdir="\"vp/parareal-ng-${j}-k-${i}-\""
    done
done