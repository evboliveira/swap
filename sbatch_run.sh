#!/bin/bash
#SBATCH -J DMRG-PlaRot-Sample
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=100G
#SBATCH --account=rrg-pnroy
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=estevao.deoliveira@uwaterloo.ca
#SBATCH --array=51-100

module load julia/1.12.5

<<<<<<< HEAD
julia /home/evbdeoli/Documents/swap/main.jl $SLURM_ARRAY_TASK_ID
=======
julia /home/evbdeoli/Documents/swap/sample.jl $SLURM_ARRAY_TASK_ID
>>>>>>> f036c41a5afc1498d72e4ba88bb84c2961366f4c
