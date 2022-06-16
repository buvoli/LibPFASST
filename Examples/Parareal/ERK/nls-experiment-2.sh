#!/bin/bash

# =============================================================================
#  Runs the second NLS experiment from paper
# =============================================================================

RUNNER='mpirun'
FLAGS_ERK='-n 1'
FLAGS_PARA='-n 1'

make # ensure that code is compiled

set -x # (print commands are they run)

# Serial ERK (Reference)
$RUNNER $FLAGS_ERK main.1d.exe params/NLS-serial-fine.nml ic_type=2 nsteps_rk=524288 rho=0 outdir='"nls-exp2/serial-ref-"'

# Serial ERK (Fine)
$RUNNER $FLAGS_ERK main.1d.exe params/NLS-serial-fine.nml ic_type=2 nsteps_rk=65536 rho=0.0245436926061703 outdir='"nls-exp2/serial-fine-"'

# Parareal ERK K=0,..,6
for i in {0..6}
do
    $RUNNER $FLAGS_PARA main.1d.exe params/NLS-parareal.nml ic_type=2 niters="${i}" rho=0.0245436926061703 outdir="\"nls-exp2/parareal-k-${i}-\""
done