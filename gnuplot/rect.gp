set term postscript eps enhanced size 8cm, 8cm color
set output "rdf_emim_dca_v2.eps" 
set encoding iso_8859_1	

set multiplot layout 3,1
set lmargin at screen 0.13
set rmargin at screen 0.98

set xrange [0:30]

#################################################################################
# plot 1
#################################################################################
set tmargin at screen 0.90
set bmargin at screen 0.65

set yrange [0:2.5]
set ylabel 'g(r)' offset 1,0

set x2tics
set xtics ( "" 0, "" 5, "" 10, "" 15, "" 20, "" 25, "" 30)
set x2label 'r / \305'
set label 1 "a) C_2mim - N(CN)_2" at graph 0.6,0.8
set label 2 "s=1" at graph 0.15,0.28 font ",size 12"
set label 3 "s=2" at graph 0.34,0.28 font ",size 12"
set label 4 "s=3" at graph 0.55,0.28 font ",size 12"
set label 5 "s=4" at graph 0.77,0.28 font ",size 12"

set key left 

plot \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_anion_g000.dat" w l lc rgb "#0063a6" lw 2.5 title "260K", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_anion_shell1_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_anion_shell2_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_anion_shell3_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_anion_shell4_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "" , \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_anion_g000.dat" w l lc rgb "#dd4814" lw 2.5 title "340K", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_anion_shell1_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_anion_shell2_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_anion_shell3_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_anion_shell4_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title ""\

#################################################################################
# plot 2
#################################################################################
set tmargin at screen 0.62
set bmargin at screen 0.38

set yrange [0:2.1]
set ytics ("0" 0, "0.5" 0.5, "1" 1, "1.5" 1.5, "2" 2)
unset x2label

set ylabel 'g(r)' offset 0,0
set x2tics ( "" 0, "" 5, "" 10, "" 15, "" 20, "" 25, "" 30)
set label 1 "b) C_2mim - C_2mim" at graph 0.6,0.8
set label 2 "s=1" at graph 0.19,0.28 font ",size 12"
set label 3 "s=2" at graph 0.38,0.28 font ",size 12"
set label 4 "s=3" at graph 0.59,0.28 font ",size 12"
set label 5 "s=4" at graph 0.81,0.28 font ",size 12"

plot \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_cation_g000.dat" w l lc rgb "#0063a6" lw 2.5 title "260 K", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_cation_shell1_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_cation_shell2_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_cation_shell3_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_cation_cation_shell4_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_cation_g000.dat" w l lc rgb "#dd4814" lw 2.5 title "340 K", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_cation_shell1_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_cation_shell2_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_cation_shell3_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_cation_cation_shell4_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title ""
\
#################################################################################
# plot 3
#################################################################################
set tmargin at screen 0.35
set bmargin at screen 0.10

set yrange [0:1.9]
set ylabel 'g(r)'
set ytics ("0" 0, "0.5" 0.5, "1" 1, "1.5" 1.5)


set xlabel 'r / \305'

set xtics ("0" 0, "5" 5, "10" 10, "15" 15, "20" 20, "25" 25, "30" 30)
set label 1 "c) N(CN)_2 - N(CN)_2" at graph 0.6,0.8
set label 2 "s=1" at graph 0.07,0.17 font ",size 12"
set label 3 "s=2" at graph 0.25,0.28 font ",size 12"
set label 4 "s=3" at graph 0.50,0.28 font ",size 12"
set label 5 "s=4" at graph 0.70,0.28 font ",size 12"
#set object 1 rect from graph 0.09, graph 0.17 to graph 0.13, graph 0.17 back
#set object 1 fc rgb 'grey' fillstyle solid 1.0
set object 1 rect from graph 0.13, graph 0.16 to graph 0.17, graph 0.17 back
set object 1 rect fc rgb "black" fillstyle solid 1.0

plot \
"/home/student1/analysis/gfunctions/rdf_dca_260K_anion_anion_g000.dat" w l lc rgb "#0063a6" lw 2.5 title "260 K", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_anion_anion_shell1_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_260K_anion_anion_shell2_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "",\
"/home/student1/analysis/gfunctions/rdf_dca_260K_anion_anion_shell3_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "",\
"/home/student1/analysis/gfunctions/rdf_dca_260K_anion_anion_shell4_g000.dat" w l lc rgb "#0063a6" lw 2.5 dashtype 3 title "",\
"/home/student1/analysis/gfunctions/rdf_dca_340K_anion_anion_g000.dat" w l lc rgb "#dd4814" lw 2.5 title "340 K", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_anion_anion_shell1_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "",\
"/home/student1/analysis/gfunctions/rdf_dca_340K_anion_anion_shell2_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "",\
"/home/student1/analysis/gfunctions/rdf_dca_340K_anion_anion_shell3_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "", \
"/home/student1/analysis/gfunctions/rdf_dca_340K_anion_anion_shell4_g000.dat" w l lc rgb "#dd4814" lw 2.5 dashtype 3 title "" 
\
