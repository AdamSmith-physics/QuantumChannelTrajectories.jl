#!/bin/bash

# If logs directory does not exist, create it
if [ ! -d "logs" ]; then
    mkdir logs
fi

dt=0.11
p=0.22
Nx=4
Ny=4
V=0.0
b=0.0
num_iterations=100
steps=100
fermions=true 

max=15
for i in `seq 1 $max`
do
    nohup julia --project=. run/run_trajectories.jl $i $dt $p $Nx $Ny $V $b $num_iterations $steps $fermions > logs/log_$i.out 2>&1 &
    echo "Process ID for run $i:" 
    echo $!
done