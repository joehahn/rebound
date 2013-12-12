#!/bin/gnuplot
set output "J2precession.pdf"
set key top left
set terminal pdf monochrome dashed enhanced size 3in,2in
set xlabel "time [years]"
set ylabel "pericenter [deg]"
set autoscale xfix 
set yrange [0:360]
mod360(x) = (x>360.)?mod360(x-360.):((x<0.)?mod360(x+360.):x)
n0 = 792.4350417074
wdot = -n0*(3./2.*16298e-6/3./3.)
set ytics 90

plot \
"a.txt" u ($1/2./pi):(mod360($4/pi*180)) w p t "simulation", \
mod360(wdot*(x*2.*pi)/pi*180.) w l t "linear theory"

