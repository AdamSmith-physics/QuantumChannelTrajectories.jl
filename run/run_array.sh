#!/bin/bash

# defq, shortq, hmemq
#SBATCH --partition=defq
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2g
#SBATCH --time=08:00:00

# Put slurm outputs in a logs directory to keep working directory clean
#SBATCH -o ./logs/output-%A_%a.out # STDOUT

#below use Linux commands, which will run on compute node. This needs to be specific to your #application

TASK=${SLURM_ARRAY_TASK_ID}

dt=0.25
p=0.25
Nx=4
Ny=4
V=1.0
b=0.1
num_iterations=500
steps=100
fermions=false 

echo "Running Job $TASK on `hostname`"
cd ${SLURM_SUBMIT_DIR}

module load julia-uoneasy/1.10.4-linux-x86_64

julia --project=. run/run_trajectories.jl $TASK $dt $p $Nx $Ny $V $b $num_iterations $steps $fermions

echo "Finished job now"