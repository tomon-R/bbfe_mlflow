set datafile separator ","
set xlabel "Time"
set ylabel "Center of mass"
set xrange [0.0 : 3.0]
set yrange [0.4 : 1.5]
set mxtics 5
set mytics 5
set key bottom left reverse
set colorsequence classic
set terminal pngcairo
set output "output.png"
plot \
"data/3d-bubble/Safi_TC1_center_LBM_h256.csv" index 0 using 1:2 w p ps 1.25 pt 4 lt 1 title "LBM (Safi, et. al., 1/h=256)", \
"results/bubble_N0/result_bubble.csv" index 0 using 1:2 w l lt 2 title "Present N1 (1/h=40)", \
"results/bubble_N2/result_bubble.csv" index 0 using 1:2 w l lt 3 title "Present N2 (1/h=40)", \
"results/bubble_N4/result_bubble.csv" index 0 using 1:2 w l lt 4 title "Present N4 (1/h=40)", \
"results/bubble_N8/result_bubble.csv" index 0 using 1:2 w l lt 5 title "Present N8 (1/h=40)", \