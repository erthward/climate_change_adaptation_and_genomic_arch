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
echo "NOW PLOTTING PHENOTYPIC SHIFT."
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_change_phenotypic_space_ch2.py > pheno_plotting.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_paneled_phenotypic_shift_image.py > pheno_panelling.pyout

# make plots for population size, average fitness, and average phenotype
echo "NOW PLOTTING Nt, mean_fit, AND mean_z OVER TIME."
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_Nt_mean_fit_mean_z_over_time.py > Nt_mean_fit_mean_z_plotting.pyout


# plot directional gene flow results
echo "NOW FITTING ALL PIDs VON MISES DISTRIBUTIONS PARAMS"
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/estimate_von_mises_mixture_params_all_files.py > gene_flow_param_fitting.pyout
echo "NOW CONCATENATING FITTED VON MISES PARAMS CSVs"
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/concatenate_fitted_von_mises_mixture_params_files.py > gene_flow_param_concat.pyout
echo "NOW PLOTTING FITTED GENE FLOW DIRECTIONAL DISTRIBUTIONS."
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/vis_vM_mix_dists.py > gene_flow_plotting.pyout

# run stats tests and output results

# copy all results to BDrive
echo "NOW PUSHING EVERYTHING TO BDRIVE USING RCLONE"
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/phenotypic_shift_*.png`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
rclone copy /global/scratch/users/drewhart/ch2/output/analysis/ch2_gene_flow_dir_analysis.png bdrive:ch2_outputs/analysis/
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/ch2_*_over_time.jpg`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
