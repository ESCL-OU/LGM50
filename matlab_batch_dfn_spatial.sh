#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=results/dfn_spatial_%j_stdout.txt
#SBATCH --error=results/dfn_spatial_%j_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=dfn_spatial
#SBATCH --mail-user=your_email@example.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --chdir=/home/yuichi/SPM-OP

module load MATLAB/2022a
mkdir -p results results/DFN
matlab -nodisplay -nosplash -nodesktop -batch "run_test_dfn_spatial_batch"
