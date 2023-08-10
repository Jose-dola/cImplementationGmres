unset key
set size square

#set xlabel 'x'; set ylabel 'y'; set zlabel 'error'
#splot 'kk' u 1:2:3 w d


#pause -1
set xrange [-1:1]; set yrange [-1:1]
set xtics -1,0.5,1
set ytics -1,0.5,1
plot 'convergen.txt' u 1:2 w p, -x, 'art' 

#pause -1
#
#
#set term postscript eps enhanced color solid "Helvetica" 24
#set out 'kk.eps'
#replot
#set term x11

