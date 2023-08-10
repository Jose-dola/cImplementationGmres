
unset key
set xlabel "m"
set ylabel "Flops"

plot for[i=0:9] 'kk0p9' u 2:3 ev :::i::i w lp

set term postscript eps enhanced color solid "Helvetica" 24
set out 'ej258jordan2.eps'
replot
