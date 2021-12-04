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

python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/plot_change_phenotypic_space_ch2.py > pheno_plot.pyout
module load rclone
for f in `ls ./phenotypic_shift_L*_G*.png`; do rclone $f bdrive:ch2_outputs; done

