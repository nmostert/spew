#!/bin/bash
X=6
for v in 5000 10000 15000 20000 25000 30000 35000 40000
do
./tephra2-2012 inputs/plume_height_trial_$X/plume_height_trial_$v.conf inputs/plume_height_trial_$X/plume_height_trial_grid.utm inputs/noWind > inputs/plume_height_trial_$X/plume_height_trial_$v.txt 
done
