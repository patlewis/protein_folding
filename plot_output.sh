#!/bin/bash

# This script only graphs the lowest-energy structure of the folded protein.
# It could be modified to graph a few of them, but I really don't want to 
# write that.  

if [ "$#" -ne 2 ]; then
    echo "usage: ./plot_output.sh <data_file> <output_plot.imagefile>"
    exit
fi


coordinates=$(grep Coordinates $1 -m 1 | awk '{print $2}')
names=$(grep Names $1 -m 1 | awk '{print $2}')
echo $coordinates
echo $names

./graphmaker.m $coordinates $2
