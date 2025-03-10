set datafile separator ","
set xlabel "t \sqrt(2g/L)"
set ylabel "Z/L"
set xrange [0.0 : 3.0]
set yrange [1.0 : 4.0]
set mxtics 5
set mytics 5
set key bottom Right reverse font ",8"
set colorsequence classic
set terminal pngcairo
set output "output.png"
plot \
"data/dambreak/martin_a_1125_n_sqrt2.csv" index 0 using ($2):($1) w p ps 1.25 pt 4 lt 2 title "Exp. (Martin & Moyce, 1.125 in.)",\
"data/dambreak/martin_a_225_n_sqrt2.csv" index 0 using ($2):($1) w p ps 1.25 pt 6 lt 3 title "Exp. (Martin & Moyce, 2.25 in.)",\
"data/dambreak/Koshizuka_exp_1995.csv" index 0 using ($1):($2) w p ps 1.25 pt 10 lt 4 title "Exp. (Koshizuka et al.)",\
"data/dambreak/Koshizuka_MPS_1996.csv" index 0 using ($1):($2) w p ps 1.25 pt 8 lt 7 title "MPS (Koshizuka et al.)",\
"data/dambreak/noslip_edge_position_a.csv" index 0 using ($1):($2) w p ps 1.25 pt 12 lt 7 title "VOF noslip a (Ota et al.)",\
"data/dambreak/noslip_edge_position_b.csv" index 0 using ($1):($2) w p ps 1.25 pt 15 lt 7 title "VOF noslip b (Ota et al.)",\
"results/dambreak_N0_40_6_40/result_dambreak.csv" index 0 using ($1 * sqrt(2*9.81/0.146)):($2/0.146) w l lt 1 title "Present (dx=0.0146)"