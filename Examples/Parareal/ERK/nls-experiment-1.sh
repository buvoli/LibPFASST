#!/bin/bash

make USE_FFT=TRUE USE_FFTW=TRUE

# # Serial ERK (Reference)
./main.1d.exe params/NLS-serial-fine.nml ic_type=1 nsteps_rk=524288 rho=0 outdir='"nls-exp1/serial-ref"'

# Serial ERK (Fine)
./main.1d.exe params/NLS-serial-fine.nml ic_type=1 nsteps_rk=65536 rho=0.0245436926061703 outdir='"nls-exp1/serial-fine"'

# Parareal ERK K=0,..,6
for i in {0..6}
do
    mpirun -n 8 --oversubscribe ./main.1d.exe params/NLS-parareal.nml ic_type=1 niters="${i}" rho=0.024543692606170 outdir="\"nls-exp1/parareal-k-${i}\""
done