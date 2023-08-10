
unset key

set xlabel "m"
set ylabel "Flops"
#beta=i*0.01;
#print i,beta
#plot 'kk' u (log10($2)):4 ev :::i::i w lp, '' u (log10($2)):6 ev :::i::i w lp, 'kk' u (log10($2)):8 ev :::i::i w lp, '' u (log10($2)):10 ev :::i::i w lp, 'kk' u (log10($2)):12 ev :::i::i w lp, '' u (log10($2)):14 ev :::i::i w lp, 'kk' u (log10($2)):16 ev :::i::i w lp, '' u (log10($2)):18 ev :::i::i w lp, 'kk' u (log10($2)):20 ev :::i::i w lp
#i=i+1; reread;


plot for[i=0:10] 'kk0p9' u 2:3 ev :::i::i w lp


set term postscript eps enhanced color solid "Helvetica" 24
set out 'ej258.eps'
replot
