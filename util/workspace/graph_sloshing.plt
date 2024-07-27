set datafile separator ","
set xlabel "Time [s]"
set ylabel "Height [m]"
set xrange [0.0 : 10]
set yrange [0.3 : 0.9]
set mxtics 5
set mytics 5
set key bottom Left reverse
set colorsequence classic
plot \
"data/sloshing/result_sloshing.csv" index 0 using ($1):($2) w p ps 1.25 pt 3 lt 1 title "Present (dx=0.025)",\
"data/sloshing/Okamoto1992.csv" index 0 using ($1):($2) w p ps 1.25 pt 4 lt 2 title "Exp. (Okamoto, 1992))",\
"data/sloshing/Sakuraba2001.csv" index 0 using ($1):($2) w p ps 1.25 pt 8 lt 7 title "Sim. (Sakuraba 2001)"