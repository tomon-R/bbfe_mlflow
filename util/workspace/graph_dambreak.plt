set datafile separator ","
set xlabel "t \sqrt(2g/L)"
set ylabel "Z/L"
set xrange [0.0 : 3.0]
set yrange [1.0 : 4.0]
set mxtics 5
set mytics 5
set key bottom Left reverse
set colorsequence classic
plot \
"data/dambreak/result_dambreak.csv" index 0 using ($1 * sqrt(2*9.81/0.146)):($2/0.146) w p ps 1.25 pt 3 lt 8 title "Present (dx=0.0146)",\
"data/dambreak/martin_a_1125_n_sqrt2.csv" index 0 using ($2):($1) w p ps 1.25 pt 4 lt 1 title "Exp. (Martin & Moyce, 1.125 in.)",\
"data/dambreak/martin_a_225_n_sqrt2.csv" index 0 using ($2):($1) w p ps 1.25 pt 6 lt 3 title "Exp. (Martin & Moyce, 2.25 in.)"