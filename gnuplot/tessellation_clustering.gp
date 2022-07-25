set term pngcairo enhanced color font 'Helvetica,10'
set encoding iso_8859_1
set out "tessellation_clustering_rgb3.png"

unset key
set yrange[-5:205]
set xlabel "time / ns"
set ylabel "Micellar size"

#set palette negative rgb 23,28,3
#set palette defined ( 0 "green", 1 "blue", 2 "red", 3 "purple")
set palette file "rgb_values3.txt"
#show palette palette 200
set pm3d map interpolate 0,0 

#splot "voronoi_clustering.dat"
splot "12k_random.dat"
