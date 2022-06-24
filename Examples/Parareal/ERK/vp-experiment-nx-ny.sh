#!/bin/bash

# =============================================================================
#  Runs the VP experiment from paper with variable nx,ny for testing
# =============================================================================

# RUNNER='mpirun'
# FLAGS_ERK='-n 1'
# FLAGS_PARA='-n 4'

RUNNER='srun'
FLAGS_ERK='--nodes=1 --ntasks=1'
FLAGS_PARA='--nodes=64 --ntasks=2048 --ntasks-per-node=32 --cpus-per-task=1'

nx=128
ny=128
base_dir="vp-nx-${nx}-ny-${ny}"

make DIM=2 # ensure that code is compiled

set -x # (print commands are they run)

# Serial ERK (Reference)
if [ ! -d "dat/${base_dir}/serial-ref-P0001/" ] ; then
    $RUNNER $FLAGS_ERK main.2d.exe params/VP-serial-fine.nml nsteps_rk=262144 nx="${nx}" ny="${ny}" rho=0 outdir="\"${base_dir}/serial-ref-\""
fi

# Serial ERK (Fine)
if [ ! -d "dat/${base_dir}/serial-fine-P0001/" ] ; then
    $RUNNER $FLAGS_ERK main.2d.exe params/VP-serial-fine.nml nsteps_rk=131072 nx="${nx}" ny="${ny}" rho=0.0245436926061703 outdir="\"${base_dir}/serial-fine-\""
fi

# Parareal ERK NG=1,...,3 and K=12
for j in {2..3}
do
    # Serial ERK (Coarse)
    nsteps_rk=$(( j * 2048 ))
    if [ ! -d "dat/${base_dir}/serial-coarse-ng-${j}-P0001/" ] ; then
        $RUNNER $FLAGS_ERK main.2d.exe params/VP-serial-coarse.nml nsteps_rk="${nsteps_rk}" nx="${nx}" ny="${ny}" rho=0.0245436926061703 outdir="\"${base_dir}/serial-coarse-ng-${j}-\""
    fi

    #Parareal K=12
    $RUNNER $FLAGS_PARA main.2d.exe params/VP-parareal.nml nx="${nx} ${nx}" ny="${ny} ${ny}" niters=12 nsteps_rk="${j} 64" rho=0.0245436926061703 outdir="\"${base_dir}/parareal-ng-${j}-\""
done