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
## Command(s) to run:

# load needed modules
module load python r gsl gcc rclone

# plot phenotypic shift results, then push results to BDrive
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_change_phenotypic_space_ch2.py > pheno_plotting.pyout
for f in `ls ./phenotypic_shift_L*_G*.png`; do rclone copy $f bdrive:ch2_outputs/analysis/; done

# plot directional gene flow results
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/estimate_von_mises_mixture_params_all_files.py > gene_flow_param_fitting.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/concatenate_fitted_von_mises_mixture_params_files.py > gene_flow_param_concat.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/vis_vM_mix_dists.py > gene_flow_plotting.pyout
rclone copy ch2_gene_flow_dir_analysis.png bdrive:ch2_outputs/analysis/

# make plots for population size and average fitness

# run stats tests and output results
