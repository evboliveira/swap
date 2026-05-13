#!/bin/bash
#SBATCH -J DMRG-PlaRot
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=100G
#SBATCH --account=rrg-pnroy
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=estevao.deoliveira@uwaterloo.ca
#SBATCH --array=2-100

module load julia/1.11.3

julia /home/evbdeoli/swap/main.jl $SLURM_ARRAY_TASK_ID