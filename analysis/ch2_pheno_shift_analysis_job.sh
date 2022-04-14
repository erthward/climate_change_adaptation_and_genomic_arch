#!/bin/bash
# Job name:
#SBATCH --job-name=ch2_anal
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

# plot phenotypic shift results
echo "NOW PLOTTING PHENOTYPIC SHIFT."
  # low redundancy
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_change_phenotypic_space.py lo non-null > pheno_shift_loREDUND.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_paneled_image.py SCAT lo non-null > pheno_panelling_loREDUND.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_change_phenotypic_space.py lo null > pheno_shift_NULL_loREDUND.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_paneled_image.py SCAT lo null > pheno_panelling_NULL_loREDUND.pyout
  # high redundancy
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_change_phenotypic_space.py hi non-null > pheno_shift_hiREDUND.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_paneled_image.py SCAT hi non-null > pheno_panelling_hiREDUND.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_change_phenotypic_space.py hi null > pheno_shift_NULL_hiREDUND.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_paneled_image.py SCAT hi null > pheno_panelling_NULL_hiREDUND.pyout

# copy all results to BDrive
echo "NOW PUSHING EVERYTHING TO BDRIVE USING RCLONE"
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/phenotypic_shift_*.jpg`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
