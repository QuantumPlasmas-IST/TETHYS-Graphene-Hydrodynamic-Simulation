set hidden3d

# axis 
set xlabel "time"
set ylabel "x"
set zlabel "density"

# plot 
splot "data/sin_yM.txt" using 1:2:3 every 3:10 ps 0.001 with lines