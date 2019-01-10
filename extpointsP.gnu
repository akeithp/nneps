set terminal postscript eps color solid "Helvetica" 16
#set terminal postscript eps color colortext
set output "./ESPP.eps"
set title "Array length vs. Search time Parallel"
set xlabel "n"
set ylabel "Time (microsecs)"
set style func linespoints
set pointsize 1
set key left
plot [] [] \
	'./RESULTSP/EPSParallel01' using 1:3 title "EPSP 1 threads" with linespoints lt 2 pt 3 lw 0.5, \
	'./RESULTSP/EPSParallel02' using 1:3 title "EPSP 2 threads" with linespoints lt 3 pt 4 lw 0.5, \
	'./RESULTSP/EPSParallel04' using 1:3 title "EPSP 4 threads" with linespoints lt 4 pt 5 lw 0.5, \
	'./RESULTSP/EPSParallel08' using 1:3 title "EPSP 8 threads" with linespoints lt 5 pt 6 lw 0.5

