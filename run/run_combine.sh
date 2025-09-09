#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1g
#SBATCH --time=0:05:00
#SBATCH --output=logs/_fer_combine_trot_%a.out

#below use Linux commands, which will run on compute node. This needs to be specific to your #application

echo "Running on `hostname`"
cd ${SLURM_SUBMIT_DIR}

module load julia-uoneasy/1.10.4-linux-x86_64

julia --project=. run/combine_data.jl $SLURM_ARRAY_TASK_ID $SLURM_JOB_ID

echo "Finished job now"