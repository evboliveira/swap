#!/bin/bash
#SBATCH -J swap_100k_samples
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --account=rrg-pnroy
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=estevao.deoliveira@uwaterloo.ca
#SBATCH --array=2-150

module load julia/1.12.5

julia /home/evbdeoli/Documents/swap/sample.jl $SLURM_ARRAY_TASK_ID