#!/bin/bash

# =============================================================================
#  Runs the VP experiment from paper
# =============================================================================

# RUNNER='mpirun'
# FLAGS_ERK='-n 1'
# FLAGS_PARA='-n 1'

RUNNER='srun'
FLAGS_ERK='--nodes=1 --ntasks=1'
FLAGS_PARA='--nodes=3 --ntasks=96 --ntasks-per-node=32 --cpus-per-task=1'

make DIM=2 # ensure that code is compiled

set -x # (print commands are they run)

# Serial ERK (Fine)
if [ ! -d "dat/vp/short-serial-fine-P0001/" ] ; then
    $RUNNER $FLAGS_ERK main.2d.exe params/VP-short-serial-fine.nml nsteps_rk=6144 rho=0.0245436926061703 outdir='"vp/short-serial-fine-"'
fi

# Parareal ERK NG=1,...,3 and K=12
for j in {2..3}
do
    # Serial ERK (Coarse)
    nsteps_rk=$(( j * 96 ))
    if [ ! -d "dat/vp/short-serial-coarse-ng-${j}-P0001/" ] ; then
        $RUNNER $FLAGS_ERK main.2d.exe params/VP-short-serial-coarse.nml nsteps_rk="${nsteps_rk}" rho=0.0245436926061703 outdir="\"vp/short-serial-coarse-ng-${j}-\""
    fi

    # Parareal ERK K=13
    $RUNNER $FLAGS_PARA main.2d.exe params/VP-short-parareal.nml predictor_proc_group_size="32" niters="28" nsteps_rk="${j} 64" rho=0.0245436926061703 outdir="\"vp/short-parareal-ng-${j}-\""
done