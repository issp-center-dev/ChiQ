#!/bin/bash

# NPROCS=1
# (number of MPI processes) must be defined elsewhere

if [ -z "$NPROCS" ]; then
    echo "NPROCS, number of MPI processes, is not defined"
    echo "Example: $ NPROCS=1 sh ./run.sh"
    exit 1
fi

echo "NPROCS = $NPROCS"

set -uex

# ---------------------------------
# DCore

echo ""
echo "##################"
echo " DCore"
echo "##################"
echo ""

ini=dmft_square.in
dcore_pre $ini 1>dcore_pre.log 2>dcore_pre.err
dcore $ini --np $NPROCS 1>dcore.log 2>dcore.err
dcore_check $ini 1>dcore_check.log 2>dcore_check.err

# ---------------------------------
# BSE

echo ""
echo "##################"
echo " BSE"
echo "##################"
echo ""

# Print version
chiq_main.py --version

# Generate q_path.dat
gen_qpath.py $ini qpath.in

# Calc chi_0(q)
dcore_chiq.py $ini --np $NPROCS 1>dcore_chiq.log 2>dcore_chiq.err

# SCL3
calc_Iq_scl.py scl_2pole.in

# Chiq post script
chiq_post.py bse.in 1>chiq_post.log 2>chiq_post.err

# Plot BSE results
(
	cd bse
	plot_chiq_path.py ../q_path.dat chi_q_scl_eigen.dat --mode='scl'
	plot_chiq_path.py ../q_path.dat I_q_scl_eigen.dat --mode='Iq'
)

echo "Finished successfully"
