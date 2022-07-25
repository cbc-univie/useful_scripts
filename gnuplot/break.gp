set term pdfcairo enhanced
set out "break2.pdf"

set multiplot
set border 1+2+4
set lmargin at screen 0.10
set rmargin at screen 0.26
set bmargin at screen 0.15
set tmargin at screen 0.95

set ylabel "Int."
unset ytics
set xtics ("9.9" 9.9, "9.7" 9.7) out
set yrange [-2e10:1e12]
set xrange [10:9.6]

set arrow 1 from 9.59, 0e10 to 9.61, -4e10 nohead
set arrow 2 from 4.09, 0e10 to 4.11, -4e10 nohead
set arrow 3 from 9.59, 101e10 to 9.61, 99e10 nohead
set arrow 4 from 4.09, 101e10 to 4.11, 99e10 nohead

plot "spec3.txt" u 1:2 w l lc rgb 'black' notitle

unset ytics
unset ylabel
set border 1+4+8
set xlabel "{/Symbol d} / ppm" offset -5,0
set yrange [-2e10:]
set lmargin at screen 0.29
set rmargin at screen 0.93
set bmargin at screen 0.15
set tmargin at screen 0.95
set xrange [4.1:2.5]
set xtics ("3.8" 3.8, "3.6" 3.6, "3.4" 3.4, "3.2" 3.2, "3.0" 3.0, "2.8" 2.8, "2.6" 2.6)
#set arrow 1 from 4.5*pi-1e-3, -1-1e-3 to 4.5*pi+1e-3, -1+1e-3 nohead
#set arrow 2 from 4.5*pi-1e-3, 1-1e-3 to 4.5*pi+1e-3, 1+1e-3 nohead
plot "spec3.txt" u 1:2 w l lc rgb 'black' t "C_4H_8O_2"
unset multiplot
