#use: ./runCrayIntera.sh 96

aprun -n $1 ./shapes #> network_output 2>&1

#cp activity_0.csv ../plot_activity/activity_0.csv
#cp activity_1.csv ../plot_activity/activity_1.csv

