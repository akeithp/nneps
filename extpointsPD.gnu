set terminal postscript eps color solid "Helvetica" 16
#set terminal postscript eps color colortext
set output "./EPSPDAC.eps"
set title "Array length vs. Search time Parallel DAC, Normal Distribution"
set xlabel "n"
set ylabel "Time (microsecs)"
set style func linespoints
set pointsize 1
set key left
plot [] [] \
	'./RESULTSP/EPSParallelD11' using 1:3 title "EPSP 1 threads" with linespoints lt 2 pt 3 lw 0.5, \
	'./RESULTSP/EPSParallelD12' using 1:3 title "EPSP 2 threads" with linespoints lt 3 pt 4 lw 0.5, \
	'./RESULTSP/EPSParallelD14' using 1:3 title "EPSP 4 threads" with linespoints lt 4 pt 5 lw 0.5, \
	'./RESULTSP/EPSParallelD18' using 1:3 title "EPSP 8 threads" with linespoints lt 5 pt 6 lw 0.5

