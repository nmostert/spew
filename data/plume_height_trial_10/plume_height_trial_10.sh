#!/bin/bash
X=10
for v in 10000 11500 13000 14500 16000 17500 19000 20500 22000 23500 25000 26500 28000 29500 31000 32500 34000 35500 37000 38500 40000
do
./tephra2-2012 inputs/plume_height_trial_$X/plume_height_trial_$v.conf inputs/plume_height_trial_$X/plume_height_trial_grid.utm inputs/noWind > inputs/plume_height_trial_$X/plume_height_trial_$v.txt 
done
