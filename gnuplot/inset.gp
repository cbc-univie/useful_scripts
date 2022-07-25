set term postscript eps enhanced size 8cm, 6cm color 
set output "taud2.eps"

set xlabel "D^{-1} / 10^{11} m^{-2} s"
set ylabel "{/Symbol t}_{sp} / ps"

set key at screen 0.43,0.9
set ytics('0' 0, '1000' 100, '2000' 200, '3000' 300, '4000' 400, '5000' 500, '6000' 600, '7000' 700, '8000' 800)
set xtics('0' 0, '2' 2e11, '4' 4e11, '6' 6e11, '8' 8e11, '10' 1e12, '12' 1.2e12)

set xrange [-8e10:1.20e12]
set yrange [-80:800]

f1(x) = k1*x + d1
f2(x) = k2*x + d2
f3(x) = k3*x + d3
f4(x) = k4*x + d4

fit [][] f1(x) "../taud_dca_cat.dat" u 1:($2*10**11) via k1, d1
fit [][] f2(x) "../taud_dca_ani.dat" u 1:($2*10**11) via k2, d2
fit [][] f3(x) "../taud_otf_cat.dat" u 1:($2*10**11) via k3, d3
fit [][] f4(x) "../taud_otf_ani.dat" u 1:($2*10**11) via k4, d4

set grid
set multiplot

plot \
f1(x) w l dt "-" lc rgb "#9ecae1" notitle, \
"../taud_dca_cat.dat" u 1:($2*10**11) ps 0.8 pt 15 lc rgb "#0063a6" t "C_2mim", \
f3(x) w l dt "-" lc rgb "#fc9772" notitle, \
"../taud_otf_cat.dat" u 1:($2*10**11) ps 0.8 pt 15 lc rgb "#dd4814" t "C_2mim", \
f2(x) w l dt "-" lc rgb "#9ecae1" notitle, \
"../taud_dca_ani.dat" u 1:($2*10**11) ps 0.8 pt 6  lc rgb "#0063a6" t "N(CN)_2", \
f4(x) w l dt "-" lc rgb "#fc9772" notitle, \
"../taud_otf_ani.dat" u 1:($2*10**11) ps 0.8 pt 8  lc rgb "#dd4814" t "OTf"

unset grid
#set lmargin at screen 0.54; set rmargin at screen 0.88
#set tmargin at screen 0.90; set bmargin at screen 0.54

set origin .50, .22
set size .4, .4
set bmargin 0; set tmargin 0; set lmargin 0; set rmargin 0
clear

unset key
unset xlabel
unset ylabel
set xrange [-8e9:1.8e11]
set yrange [-10:100]
unset ytics
#set ytics('' 0, '' 40, '' 80)
set y2tics('' 0, '300' 30, '600' 60, '' 90) offset -5 mirror
set x2tics('0' 0, '0.5' 5e10, '1.0' 1e11, '1.5' 1.5e11) offset 0,-2 mirror

unset xtics

plot \
f1(x) w l dt "-" lc rgb "#9ecae1" notitle, \
"../taud_dca_cat.dat" u 1:($2*10**11) ps 0.8 pt 15 lc rgb "#0063a6" t "C_2mim", \
f3(x) w l dt "-" lc rgb "#fc9772" notitle, \
"../taud_otf_cat.dat" u 1:($2*10**11) ps 0.8 pt 15 lc rgb "#dd4814" t "C_2mim", \
f2(x) w l dt "-" lc rgb "#9ecae1" notitle, \
"../taud_dca_ani.dat" u 1:($2*10**11) ps 0.8 pt 6  lc rgb "#0063a6" t "N(CN)_2", \
f4(x) w l dt "-" lc rgb "#fc9772" notitle, \
"../taud_otf_ani.dat" u 1:($2*10**11) ps 0.8 pt 8  lc rgb "#dd4814" t "OTf"

unset multiplot
