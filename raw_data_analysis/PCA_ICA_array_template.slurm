#!/bin/bash
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -n 1
#SBATCH -t 180          # Runtime in minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=25000           # Memory
#SBATCH -o "datapath/pca_ica_%%A_%%a.out"  # File to which STDOUT will be written, %%j inserts jobid
#SBATCH -e "datapath/pca_ica_%%A_%%a.err"  # File to which STDERR will be written, %%j inserts jobid
#SBATCH --array=1-%d

module load matlab/R2018a-fasrc01
matlab -sd "~/" -nosplash -nodesktop -r "cd 'datapath'; ICA_PCA_array($SLURM_ARRAY_TASK_ID)"