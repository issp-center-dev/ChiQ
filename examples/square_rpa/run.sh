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
bse_tool.py --version

# Generate q_path.dat
gen_qpath.py $ini qpath.in

# Calc chi_0(q)
dcore_bse $ini --np $NPROCS 1>dcore_bse.log 2>dcore_bse.err

# BSE main script
mpirun -np $NPROCS bse_tool.py bse.in 1>bse_tool.log 2>bse_tool.err

# BSE post script
bse_post.py bse.in 1>bse_post.log 2>bse_post.err

# Plot BSE results
(
	cd bse
	plot_chiq_path.py ../q_path.dat chi0_q_eigen.dat --mode='chi0'
	plot_chiq_path.py ../q_path.dat chi_q_rpa_eigen.dat --mode='rpa'
)
