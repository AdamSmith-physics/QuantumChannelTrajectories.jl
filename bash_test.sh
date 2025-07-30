#!/bin/bash

# If logs directory does not exist, create it
if [ ! -d "logs" ]; then
    mkdir logs
fi

max=10
for i in `seq 1 $max`
do
    nohup julia --project=. run/run_trajectories.jl $i > logs/log_$i.out &
done