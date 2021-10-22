#!/bin/bash
# Job name:
#SBATCH --job-name=ch2_data
#
# Account:
#SBATCH --account=fc_landgen
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
# Set amount of memory per node (since OOM errors occurred)
#SBATCH --mem-per-cpu=2500
#
# Email options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drew.hart@berkeley.edu
#
## Command(s) to run (example):

module load python gsl gcc

python /global/home/users/drewhart/genarch_and_envchange/climate_change_adaptation_and_genomic_arch/parallelize_batches.py > /global/scratch/users/drewhart/ch2/output/10-21-21/ch2_job.pyout

