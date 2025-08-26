#!/bin/bash

# defq, shortq, hmemq
#SBATCH --partition=defq
#SBATCH --array=1-18
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4g
#SBATCH --time=02:00:00
#SBATCH --output=logs/out_test_run_%a.out



SCRIPT=run/run_trajectories.jl


echo "Running Job $TASK on `hostname`"
cd ${SLURM_SUBMIT_DIR}
module load julia-uoneasy/1.10.4-linux-x86_64

julia --project=. $SCRIPT $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID

#### to run with user installed julia app
## srun ~/julia-1.10.1/bin/julia --project=. --threads 1 $SCRIPT $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID

echo "Finished job now"

