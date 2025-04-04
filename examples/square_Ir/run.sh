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

# Generate q_fbz.dat
gen_allq.py $ini

# Calc chi_0(q) and vertex part (heavy)
dcore_chiq.py $ini --np $NPROCS 1>dcore_chiq.log 2>dcore_chiq.err

# BSE main script
mpirun -np $NPROCS chiq_main.py bse.in 1>chiq_main.log 2>chiq_main.err

# Calc I(q)
calc_Iq.py  --remove 1

# I(q) -> I(r)
chiq_fft.py --input_dname I_q --output_dname I_r $ini

# BSE post script
chiq_post.py bse.in 1>chiq_post.log 2>chiq_post.err

# Plot I(r)
(
	cd bse
	plot_Ir.py --label ../label2.in I_r_eigen.dat
)

# Generate q_path.dat (only for plot)
gen_qpath.py $ini qpath.in

# Plot chi(q), I(q)
(
	cd bse
	plot_chiq_path.py ../q_path.dat --label ../label1.in chi_q_eigen.dat
	plot_chiq_path.py ../q_path.dat --label ../label2.in I_q_eigen.dat --mode='Iq'
)

echo "Finished successfully"
