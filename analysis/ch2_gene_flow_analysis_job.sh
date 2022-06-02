#!/bin/bash
# Job name:
#SBATCH --job-name=ch2A_gf
#
# Account:
#SBATCH --account=fc_landgen
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Email options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drew.hart@berkeley.edu
#
## Command(s) to run:

# load needed modules
module load python r gsl gcc rclone

# plot directional gene flow results
echo "NOW FITTING ALL PIDs VON MISES DISTRIBUTIONS PARAMS"
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/estimate_von_mises_mixture_params_all_files.py > gene_flow_param_fitting.pyout
echo "NOW CONCATENATING FITTED VON MISES PARAMS CSVs"
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/concatenate_fitted_von_mises_mixture_params_files.py > gene_flow_param_concat.pyout
echo "NOW PLOTTING FITTED GENE FLOW DIRECTIONAL DISTRIBUTIONS."
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/viz_vM_mix_dists.py > gene_flow_plotting.pyout

# copy all results to BDrive
echo "NOW PUSHING EVERYTHING TO BDRIVE USING RCLONE"
rclone copy /global/scratch/users/drewhart/ch2/output/analysis/ --include  "ch2_gene_flow_dir_analysis*" bdrive:ch2_outputs/analysis
