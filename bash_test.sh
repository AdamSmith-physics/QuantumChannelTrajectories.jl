#!/bin/bash

# If logs directory does not exist, create it
if [ ! -d "logs" ]; then
    mkdir logs
fi

dt=0.25
p=0.5
Nx=4
Ny=4
V=1.5
b=0.5
num_iterations=500
steps=200
fermions=false 

max=15
for i in `seq 1 $max`
do
    nohup julia --project=. run/run_trajectories.jl $i $dt $p $Nx $Ny $V $b $num_iterations $steps $fermions > logs/log_$i.out &
done