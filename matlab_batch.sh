#!/bin/bash
# Andrew H. Fagg
#
# Example with an array of experiments
#  The --array line says that we will execute 4 experiments (numbered 0,1,2,3).
#   You can specify ranges or comma-separated lists on this line
#  For each experiment, the SLURM_ARRAY_TASK_ID will be set to the experiment number
#   In this case, this ID is used to set the name of the stdout/stderr file names
#   and is passed as an argument to the python program
#
#
# When you use this batch file:
#  Change the email address to yours! (I don't want email about your experiments)
#  Change the chdir line to match the location of where your code is located
#
# Reasonable partitions: debug_5min, debug_30min, normal, disc_dual_a100_students
#

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
# The %j is translated into the job number
#SBATCH --output=results/spm_op_%j_stdout.txt
#SBATCH --error=results/spm_op_%j_stderr.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=spm_op_multi
#SBATCH --mail-user=yuichi@ou.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --chdir=/home/yuichi/SPM-OP
#
#################################################

module load MATLAB/2022a
mkdir -p results
matlab -nodisplay -nosplash -nodesktop -batch "run_test_spme_batch"

