unset key

plot for[i=0:0] 'kk' u 2:3 ev :::i::i w lp

set term postscript eps enhanced color solid "Helvetica" 24
set out 'krein3.eps'
replot
