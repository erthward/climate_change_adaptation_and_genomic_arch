#!/bin/bash
# Job name:
#SBATCH --job-name=ch2A_Ntfit
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
# Email options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drew.hart@berkeley.edu
#
## Command(s) to run:

# load needed modules
module load python r gsl gcc rclone

# make plots for population size, average fitness, and average phenotype
echo "NOW PLOTTING Nt, mean_fit, AND mean_z OVER TIME."
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_Nt_mean_fit_mean_z_over_time.py > Nt_mean_fit_mean_z_plotting.pyout

# make delta_Nt and delta_fit boxplots
echo "NOW PLOTTING delta_Nt and delta_fit boxplots."
Rscript --vanilla /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_boxplots.R > boxplots_and_stats.Rout

# copy all results to BDrive
echo "NOW PUSHING EVERYTHING TO BDRIVE USING RCLONE"
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/ch2_*_over_time.jpg`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/boxplot*.jpg`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
