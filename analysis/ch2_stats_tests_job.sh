#!/bin/bash
# Job name:
#SBATCH --job-name=ch2A_stat_dens
#
# Account:
#SBATCH --account=fc_landgen
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=2:00:00
#
# Email options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drew.hart@berkeley.edu
#
## Command(s) to run:

# load needed modules
module load r r-packages gsl gcc rclone

# run stats tests and output results
echo "NOW RUNNING REMAINING STATISTICAL TESTS"
Rscript --vanilla /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/run_all_stats_tests.R > test_results.Rout

# copy all results to BDrive
echo "NOW PUSHING EVERYTHING TO BDRIVE USING RCLONE"
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/test_results.Rout`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
