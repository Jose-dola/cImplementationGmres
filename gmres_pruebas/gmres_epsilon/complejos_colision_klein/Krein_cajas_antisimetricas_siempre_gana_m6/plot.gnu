unset key

set xlabel "m"
set ylabel "Flops"

plot for[i=0:10] 'kkkk' u 1:2 ev :::i::i w lp
#plot 'kkk' w lp

set term postscript eps enhanced color solid "Helvetica" 24
set out 'krein2.eps'
replot
