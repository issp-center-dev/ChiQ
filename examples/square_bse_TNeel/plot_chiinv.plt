set key left
set xlabel "T"
set ylabel "1/chi(M)"

f(x) = a*(x-TN)

fit f(x) 'chi.dat' u 1:(1.0/$2) via a, TN

plot 'chi.dat' u 1:(1.0/$2) w p pt 4 t"1/chi(M)",\
 f(x) t"a(T-T_N)",\
  0.0 lc "black" dt (10,10) t""

set term push
set size 0.6
set term postscript eps enhanced color solid
set output "chiinv.eps"
rep
set output
set size 1.0
set term pop
