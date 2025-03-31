#!/bin/bash

if [ -z "$NPROCS" ]; then
  echo "NPROCS (number of MPI processes) is not defined"
  echo "Example: $ NPROCS=1 sh ./run.sh"
  exit 1
fi

echo "NPROCS = $NPROCS"

set -uex

ROOTDIR=$(pwd)
ini=dmft_square.in

# (1) Perform preprocess
rm -f square.h5
dcore_pre $ini 1>dcore_pre.log 2>dcore_pre.err
gen_qpath.py $ini qpath.in
rm -f chi.dat
echo "# T chi" > chi.dat

# (2) Get the momentum of M point as Q_M (e.g., 16.16.00)
# Tips: In AWK, $(NF) means the last field of the line.
Q_M=$(awk '$(NF)=="M" {print $2}' < q_path.dat)

for T in 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
  echo "T = $T"
  cd "$ROOTDIR"

  # (3) Create a directory for each temperature, copy necessary files, and add T to the input file
  # NOTE: remove the directory and its contents if it already exists
  rm -rf T_$T
  mkdir -p T_$T
  for f in dmft_square.in square.h5 square_hr.dat q_path.dat bse.in; do
    cp $f T_$T
  done
  echo "T = ${T}" >> T_$T/$ini
  cd T_$T

  # (4) calculate susceptibilities by BSE
  dcore $ini --np $NPROCS 1>dcore.log 2>dcore.err
  dcore_check $ini 1>dcore_check.log 2>dcore_check.err
  dcore_chiq.py $ini --np $NPROCS 1>dcore_chiq.log 2>dcore_chiq.err
  mpirun -np $NPROCS chiq_main.py bse.in 1>chiq_main.log 2>chiq_main.err
  chiq_post.py bse.in 1>chiq_post.log 2>chiq_post.err

  # (5) get chi at Q_M and save to chi.dat
  chi=$(awk -v Q="${Q_M}" '$2==Q {print $3}' < bse/chi_q_eigen.dat)
  echo "$T $chi" >> ../chi.dat

  # Plot BSE results
  (
    cd bse
    plot_chiq_path.py ../q_path.dat chi_q_eigen.dat
  )

  cd "$ROOTDIR"
done

# Plot the inverse susceptibility
gnuplot plot_chiinv.plt

echo "Finished successfully"
