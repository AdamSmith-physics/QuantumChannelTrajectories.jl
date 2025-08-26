#!/bin/bash

# defq, shortq, hmemq
#SBATCH --partition=defq
#SBATCH --array=1-40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2g
#SBATCH --time=08:00:00
#SBATCH --output=logs/out_%j_%a.out



SCRIPT=run/run_trajectories.jl


echo "Running Job $TASK on `hostname`"
cd ${SLURM_SUBMIT_DIR}
module load julia-uoneasy/1.10.4-linux-x86_64

julia --project=. $SCRIPT $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID

#### to run with user installed julia app
## srun ~/julia-1.10.1/bin/julia --project=. --threads 1 $SCRIPT $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID

echo "Finished job now"

