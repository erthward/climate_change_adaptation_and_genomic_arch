#!/bin/bash
# Job name:
#SBATCH --job-name=ch2_sim
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

module load gsl gcc python/3.7

python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/sim/parallelize_batches.py > ch2_job.pyout

