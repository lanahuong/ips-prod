# We are loading a csv file os we specify "," as the separator
set datafile separator comma

# Load the colormap
load 'scripts/magma.pal'

# On the x-axis we display z that ranges from -20 to 20 fm
set xrange [-20 :20 ]
# On the y-axis we display r that ranges from -10 to 10 fm
set yrange [-10 : 10]
# This gives the same scale on both axis
set size ratio 0.5

# Save the graph in a file
set terminal png size 600,300
set output "tmp/density-r-z.png"

# Add title and labels
set title "Intensité de la densité nucléaire en fonction de z et r"
set xlabel "z (fm)"
set ylabel "r (fm)"
set tics nomirror out scale 1

# Plot the density matrix
plot "tmp/density-r-z.csv" matrix u (40*($1-32)/63):(20*($2-16)/31):3 with image notitle