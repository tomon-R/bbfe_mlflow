set datafile separator ","
set xlabel "Time"
set ylabel "Rise velocity"
set xrange [0.0 : 3.0]
set yrange [0.95 : 1.0]
set mxtics 5
set mytics 5
set key bottom Left reverse
set colorsequence classic
plot \
"data/3d-bubble/result_bubble.csv" index 0 using 1:4 w p ps 1.25 pt 3 lt 8 title "Present (dx=0.0146)",\
"data/3d-bubble/Safi_TC1_sphericity_LBM_h256.csv" index 0 using 1:2 w p ps 1.25 pt 4 lt 1 title "LBM (Safi, et. al., h=1/256)"