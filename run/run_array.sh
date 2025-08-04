#!/bin/bash

#SBATCH --partition=hmemq
#SBATCH --array=1-75
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40g
#SBATCH --time=3-00:00:00

# Put slurm outputs in a logs directory to keep working directory clean
#SBATCH -o ./logs/output-%A_%a.out # STDOUT

#below use Linux commands, which will run on compute node. This needs to be specific to your #application

TASK=${SLURM_ARRAY_TASK_ID}

echo "Running Job $TASK on `hostname`"
cd ${SLURM_SUBMIT_DIR}

module load julia-uoneasy/1.10.4-linux-x86_64

nohup julia --project=. run/run_trajectories.jl $TASK

echo "Finished job now"