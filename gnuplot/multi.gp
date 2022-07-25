set term pdfcairo color font ", 9"
set output "spec_all_v2.eps"

set style line 1 lc rgb '#000063a6' lw 2
set style line 2 lc rgb '#4d0063a6' lw 2
set style line 3 lc rgb '#990063a6' lw 2
set style line 4 lc rgb '#00dd4814' lw 2
set style line 5 lc rgb '#4ddd4814' lw 2
set style line 6 lc rgb '#99dd4814' lw 2
set style line 7 lc rgb '#FF8BA8' lw 2
set style line 8 lc rgb '#A3274D' lw 2
set style line 9 lc rgb '#193B63' lw 2
set style line 10 lc rgb '#7096BA' lw 2


##########################################################################
#Plot1
##########################################################################

set macros

#set xrange [0:100]
#set yrange [0:20000]


set xrange [0.00001:50]
set logscale x
set yrange [0:5]

set xlabel "{/Symbol w} / THz"
#set ylabel "{/Symbol e}^{\"}({/Symbol w})"

L2 = "set lmargin at screen 0.10; set rmargin at screen 0.50"
L1 = "set lmargin at screen 0.55; set rmargin at screen 0.95"
T1 = "set tmargin at screen 0.10; set bmargin at screen 0.27"
T2 = "set tmargin at screen 0.27; set bmargin at screen 0.44"
T3 = "set tmargin at screen 0.44; set bmargin at screen 0.61"
T4 = "set tmargin at screen 0.61; set bmargin at screen 0.78"
T5 = "set tmargin at screen 0.78; set bmargin at screen 0.95"


XL = 'set xlabel "{/Symbol w} / THz"'
NX = "unset xlabel"
XT = "set xtics ('' 0.00001, '0.0001' 0.0001, '0.001' 0.001,  '0.01' 0.01, '0.1' 0.1, '1.0' 1.0, '10' 10 )"
NT = "set xtics ('' 0.00001, '' 0.0001, '' 0.001, '' 0.01, '' 0.1, '' 1.0, '' 10 )"
NY = "unset ylabel"
YL1 = 'set ylabel "{/Symbol S}@^{\"}_{0 [260K]}"'
YL2 = 'set ylabel "{/Symbol S}@^{\"}_{0 [280K]}"'
YL3 = 'set ylabel "{/Symbol S}@^{\"}_{0 [300K]}"'
YL4 = 'set ylabel "{/Symbol S}@^{\"}_{0 [320K]}"'
YL5 = 'set ylabel "{/Symbol S}@^{\"}_{0 [340K]}"'



YR11 = "set yrange [0:6]"
YR12 = "set yrange [0:6]"
YR13 = "set yrange [0:6]"
YR14 = "set yrange [0:6]"
YR15 = "set yrange [0:6]"



YT1 = "set ytics ('' 0, '2' 2, '4' 4)"
YT2 = "set ytics ('' 0, '2' 2, '4' 4)"
YT3 = "set ytics ('' 0, '2' 2, '4' 4)"                   
YT4 = "set ytics ('' 0, '2' 2, '4' 4)"                  
YT5 = "set ytics ('' 0, '2' 2, '4' 4)"



set multiplot layout 5,2 rowsfirst

unset key

    @L1;@T1;@YT1;@YR11;@XT;@NY;
    plot \
         "/site/raid1/student2/flora/260K_trif/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 4 title "{/Symbol e}@_0^{\"}({/Symbol w})",\
         "/site/raid1/student2/flora/260K_trif/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 5 dt 4 title "{/Symbol J}@_0^{\"}({/Symbol w})",\
         "/site/raid1/student2/flora/260K_trif/analysis/new_fitprogram/gendicon_fit.dat" u 1:3 w l ls 6 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})", \
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e

@L2;@T1;@YT1;@YR11;@XT;@YL1;
plot \
     "/site/raid1/student2/flora/260K/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 1 title "{/Symbol e}^{\"}({/Symbol w})",\
     "/site/raid1/student2/flora/260K/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 2 dt 4 title "{/Symbol J}^{\"}({/Symbol w})",\
     "/site/raid1/student2/flora/260K/analysis/new_fitprogram/gendicon_fit.dat" u 1:3 w l ls 3 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})", \
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e



    @L1;@T2;@NX;@NT;@YT2;@YR12;@NY;
    plot \
"/site/raid1/student2/flora/280K_trif/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 4 title "{/Symbol e}^{\"}({/Symbol w})",\
 "/site/raid1/student2/flora/280K_trif/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 5 dt 4 title "{/Symbol J}^{\"}({/Symbol w})",\
"/site/raid1/student2/flora/280K_trif/analysis/new_fitprogram/gendicon_fit.dat" u 1:(130*$3) w l ls 6 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})",\
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e


    @L2;@T2;@NX;@NT;@YT2;@YR12;@YL2;
plot \
"/site/raid1/student2/flora/280K/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 1 title "{/Symbol e}^{\"}({/Symbol w})",\
 "/site/raid1/student2/flora/280K/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 2 dt 4 title "{/Symbol J}^{\"}({/Symbol w})",\
"/site/raid1/student2/flora/280K/analysis/new_fitprogram/gendicon_fit.dat" u 1:3 w l ls 3 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})", \
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e




   @L1;@T3;@NX;@NT;@YT3;@YR13;@NY;
    plot \
"/site/raid1/student2/flora/300K_trif/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 4 title "{/Symbol e}^{\"}({/Symbol w})",\
 "/site/raid1/student2/flora/300K_trif/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 5 dt 4 title "{/Symbol J}^{\"}({/Symbol w})",\
"/site/raid1/student2/flora/300K_trif/analysis/new_fitprogram/gendicon_fit.dat" u 1:3 w l ls 6 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})", \
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e


    @L2;@T3;@NX;@NT;@YT3;@YR13;@YL3;
plot \
"/site/raid1/student2/flora/300K/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 1 title "{/Symbol e}^{\"}({/Symbol w})",\
 "/site/raid1/student2/flora/300K/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 2 dt 4 title "{/Symbol J}^{\"}({/Symbol w})",\
"/site/raid1/student2/flora/300K/analysis/new_fitprogram/gendicon_fit.dat" u 1:3 w l ls 3 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})", \
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e



    @L1;@T4;@NX;@NT;@YT4;@YR14;@NY;
    plot \
"/site/raid1/student2/flora/320K_trif/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 4 title "{/Symbol e}^{\"}({/Symbol w})",\
 "/site/raid1/student2/flora/320K_trif/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 5 dt 4 title "{/Symbol J}^{\"}({/Symbol w})",\
"/site/raid1/student2/flora/320K_trif/analysis/new_fitprogram/gendicon_fit.dat" u 1:3 w l ls 6 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})", \
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e


    @L2;@T4;@NX;@NT;@YT4;@YR14;@YL4;
plot \
"/site/raid1/student2/flora/320K/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 1 title "{/Symbol e}^{\"}({/Symbol w})",\
 "/site/raid1/student2/flora/320K/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 2 dt 4 title "{/Symbol J}^{\"}({/Symbol w})",\
"/site/raid1/student2/flora/320K/analysis/new_fitprogram/gendicon_fit.dat" u 1:3 w l ls 3 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})", \
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e


set key top left
set key spacing 3


    @L1;@T5;@NX;@NT;@YT5;@YR15;@NY;
    set x2label "[C_2mim]OTf"
plot \
"/site/raid1/student2/flora/340K_trif/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 4 title "{/Symbol e}^{\"}({/Symbol w})",\
 "/site/raid1/student2/flora/340K_trif/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 5 dt 4 title "{/Symbol J}@_0^{\"}({/Symbol w})",\
"/site/raid1/student2/flora/340K_trif/analysis/new_fitprogram/gendicon_fit.dat" u 1:3 w l ls 6 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})", \
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e


    @L2;@T5;@NX;@NT;@YT5;@YR15;@YL5;
unset x2label; set x2label "[C_2mim]N(CN)_2"
plot \
"/site/raid1/student2/flora/340K/analysis/new_fitprogram/epsilon_fit.dat" u 1:3 w l ls 1 title "{/Symbol e}^{\"}({/Symbol w})",\
 "/site/raid1/student2/flora/340K/analysis/new_fitprogram/theta_fit.dat" u 1:3 w l ls 2 dt 4 title "{/Symbol J}@_0^{\"}({/Symbol w})",\
"/site/raid1/student2/flora/340K/analysis/new_fitprogram/gendicon_fit.dat" u 1:3 w l ls 3 dt 1 title "{/Symbol S}@_0^{\"}({/Symbol w})", \
         "-" w i lc rgb '#7a7a7a' lw 3 dt 2 notitle
         3 6
         e

@L2
set tmargin at screen 0.10; set bmargin at screen 0.95
unset xtics; unset ytics
unset x2tics; unset y2tics
unset xlabel; unset x2label
unset ylabel; unset y2label
unset key

unset logscale x

set xrange [0:7]
set yrange [0:100]

p 7*x**2.3-18 w l ls 1 dt 2

@L1
p 7*x**2.3-18 w l ls 4 dt 2

unset multiplot

