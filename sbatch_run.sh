#!/bin/bash
#SBATCH -J DMRG-PlaRot-Sample
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --account=rrg-pnroy
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=estevao.deoliveira@uwaterloo.ca
#SBATCH --array=4-150

module load julia/1.12.5

julia /home/evbdeoli/Documents/swap/sample.jl $SLURM_ARRAY_TASK_ID