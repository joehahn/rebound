#!/usr/bin/gnuplot
set output "cdf_stochastic.pdf"
set terminal pdf enhanced  dashed size 3.5in, 2.5in 

set ylabel "number of pairs"
set xlabel "period ratio"
set xrange [1:9]
set logscale x
set yrange [0:*]
set grid x
set xtics (1, 1.33, 1.5, 2.,3.,4.,5.,6.,7.,8.,9.)
set cblabel "total planet mass m_1 + m_2 [M_{jup}]"
set key top left
set cbrange [0.:.1]
plot "cdf.kepler" u 1:0 w l t "Observed KOI planets", \
"cdf.integrated.out_stochastic_1e-6" u 1:0 w l t "Migration {/Symbol t}_a=10^{4} years, {/Symbol a}=10^{-6}"  ls 4 lc -1, \
"cdf.integrated.out_stochastic_1e-5" u 1:0 w l t "Migration {/Symbol t}_a=10^{4} years, {/Symbol a}=10^{-5}"  ls 2 lc -1
