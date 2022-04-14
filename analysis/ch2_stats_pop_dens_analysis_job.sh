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
#SBATCH --time=72:00:00
#
# Email options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drew.hart@berkeley.edu
#
## Command(s) to run:

# load needed modules
module load python r gsl gcc rclone

# run stats tests and output results
echo "NOW RUNNING REMAINING STATISTICAL TESTS"
Rscript --vanilla /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/run_tests.R > test_results.Rout

# plot pop-density shift results
echo "NOW PLOTTING POP-DENSITY SHIFT."
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_change_pop_density.py > pop_dens.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_paneled_image.py DENS lo > dens_panelling_loREDUND.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_paneled_image.py DENS hi > dens_panelling_hiREDUND.pyout

# copy all results to BDrive
echo "NOW PUSHING EVERYTHING TO BDRIVE USING RCLONE"
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/pop_density_shift_*.jpg`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
