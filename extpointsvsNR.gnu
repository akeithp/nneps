set terminal postscript eps color solid "Helvetica" 16
#set terminal postscript eps color colortext
set output "./EPSvsEPSPR.eps"
set title "Array length vs. Search time, Parallel Reduction, Normal Distribution"
set xlabel "n"
set ylabel "Time (microsecs)"
set style func linespoints
set pointsize 1
set key left
plot [] [] \
	'./RESULTS/EPS1' using 1:2 title "EPS" with linespoints lt 1 pt 2 lw 0.5, \
	'./RESULTSP/EPSParallelR11' using 1:2 title "EPSPRED 1 threads" with linespoints lt 2 pt 3 lw 0.5, \
	'./RESULTSP/EPSParallelR12' using 1:2 title "EPSPRED 2 threads" with linespoints lt 3 pt 4 lw 0.5, \
	'./RESULTSP/EPSParallelR14' using 1:2 title "EPSPRED 4 threads" with linespoints lt 4 pt 5 lw 0.5, \
	'./RESULTSP/EPSParallelR18' using 1:2 title "EPSPRED 8 threads" with linespoints lt 5 pt 6 lw 0.5

