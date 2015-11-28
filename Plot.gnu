#calling this gnu file should be done as follwing
#gnuplot -e "amount = ..." Plot.gnu

set term png
set yrange[-1 to 1]
unset key
set output 'magnetisation.png'
set title 'Magnetisation of the ising system'
plot 'energyandmag.md' using 1:3 with lines

set xrange[0 to 1]
set yrange[0 to 1]
set cbrange[-1.1 to 1.1]
unset colorbox
unset key
set title "Spin Ising system"

set term png

#gnuplots does 'tot en met' 
do for [i=1:amount]{
	outfile = sprintf('Plot%03i.png',i)
	set output outfile
	plot 'Results.md' index i using 1:2:3 with points palette
}
