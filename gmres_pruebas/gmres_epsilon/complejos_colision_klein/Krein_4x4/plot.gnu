unset key
set xtics 1
set xlabel "m"
set ylabel "Flops"

#plot 'kk' u 2:3 w lp
plot for[i=0:3] 'kk' u 2:3 ev :::i::i w lp

replot
set term postscript eps enhanced color solid "Helvetica" 24
set out 'colisionkrein.eps'
replot
