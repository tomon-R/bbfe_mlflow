set datafile separator ","
set xlabel "Time"
set ylabel "Center of mass"
set xrange [0.0 : 3.0]
set yrange [0.4 : 1.5]
set mxtics 5
set mytics 5
set key bottom Left reverse
set colorsequence classic
plot \
"data/3d-bubble/result_bubble_tc2.csv" index 0 using 1:2 w l ps 1.25 pt 3 lt 8 title "Present (1/h=40)",\
"data/3d-bubble/Safi_TC2_center_LBM_h256.csv" index 0 using 1:2 w p ps 1.25 pt 4 lt 1 title "LBM (Safi, et. al., 1/h=256)",\
"data/3d-bubble/Safi_TC2_center_FeatFlow.csv" index 0 using 1:2 w p ps 1.25 pt 6 lt 1 title "FeatFlow (Safi, et. al.)"