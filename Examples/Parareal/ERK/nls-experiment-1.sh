#!/bin/bash

# =============================================================================
#  Runs the first NLS experiment from paper
# =============================================================================

RUNNER='mpirun'
FLAGS_ERK='-n 1'
FLAGS_PARA='-n 1'

make DIM=1 # ensure that code is compiled

set -x # (print commands are they run)

# Serial ERK (Reference)
if [ ! -d "dat/nls-exp1/serial-ref-P0001/" ] ; then
    $RUNNER $FLAGS_ERK main.1d.exe params/NLS-serial-fine.nml ic_type=1 nsteps_rk=524288 rho=0 outdir='"nls-exp1/serial-ref-"'
fi

# Serial ERK (Fine)
if [ ! -d "dat/nls-exp1/serial-fine-P0001/" ] ; then
    $RUNNER $FLAGS_ERK main.1d.exe params/NLS-serial-fine.nml ic_type=1 nsteps_rk=65536 rho=0.0245436926061703 outdir='"nls-exp1/serial-fine-"'
fi

# Parareal ERK K=0,..,6
for i in {0..6}
do
    $RUNNER $FLAGS_PARA main.1d.exe params/NLS-parareal.nml ic_type=1 niters="${i}" rho=0.0245436926061703 outdir="\"nls-exp1/parareal-k-${i}-\""
done