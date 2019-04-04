#!/bin/bash
X=3
for v in 0 288 3600 7200 10800 14400 18000 10000000
do
./tephra2-2012 inputs/fall_time_threshold_trial_$X/fall_time_threshold_trial_$v.conf inputs/fall_time_threshold_trial_$X/fall_time_threshold_trial_grid.utm inputs/noWind > inputs/fall_time_threshold_trial_$X/fall_time_threshold_trial_$v.txt 
done
