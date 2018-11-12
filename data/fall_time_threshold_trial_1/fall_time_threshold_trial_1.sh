#!/bin/bash
X=1
for v in 0 60 288 3600 10800 86400 259200 10000000
do
./tephra2-2012 inputs/fall_time_threshold_trial_$X/fall_time_threshold_trial_$v.conf inputs/fall_time_threshold_trial_$X/fall_time_threshold_trial_grid.utm inputs/noWind > inputs/fall_time_threshold_trial_$X/fall_time_threshold_trial_$v.txt 
done
