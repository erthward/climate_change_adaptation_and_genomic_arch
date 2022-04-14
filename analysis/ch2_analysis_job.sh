#!/bin/bash
# Job name:
#SBATCH --job-name=ch2A_all
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

# plot pop-density shift results
echo "NOW PLOTTING POP-DENSITY SHIFT."
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_change_pop_density.py > pop_dens.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_paneled_image.py DENS lo > dens_panelling_loREDUND.pyout
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_paneled_image.py DENS hi > dens_panelling_hiREDUND.pyout

# make plots for population size, average fitness, and average phenotype
echo "NOW PLOTTING Nt, mean_fit, AND mean_z OVER TIME."
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/plot_Nt_mean_fit_mean_z_over_time.py > Nt_mean_fit_mean_z_plotting.pyout

# make delta_Nt and delta_fit boxplots
echo "NOW PLOTTING delta_Nt and delta_fit boxplots."
Rscript --vanilla /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/make_boxplots.R > boxplots_and_stats.Rout

# plot directional gene flow results
echo "NOW FITTING ALL PIDs VON MISES DISTRIBUTIONS PARAMS"
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/estimate_von_mises_mixture_params_all_files.py > gene_flow_param_fitting.pyout
echo "NOW CONCATENATING FITTED VON MISES PARAMS CSVs"
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/concatenate_fitted_von_mises_mixture_params_files.py > gene_flow_param_concat.pyout
echo "NOW PLOTTING FITTED GENE FLOW DIRECTIONAL DISTRIBUTIONS."
python /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/viz_vM_mix_dists.py > gene_flow_plotting.pyout

# run stats tests and output results
echo "NOW RUNNING REMAINING STATISTICAL TESTS"
Rscript --vanilla /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/run_tests.R > test_results.Rout


# copy all results to BDrive
echo "NOW PUSHING EVERYTHING TO BDRIVE USING RCLONE"
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/pop_density_shift_*.jpg`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/phenotypic_shift_*.jpg`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/ch2_*_over_time.jpg`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
for f in `ls /global/scratch/users/drewhart/ch2/output/analysis/boxplot*.jpg`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
rclone copy /global/scratch/users/drewhart/ch2/output/analysis/ch2_gene_flow_dir_analysis*.png bdrive:ch2_outputs/analysis/
#for f in `ls /global/scratch/users/drewhart/ch2/output/output/fig_time_*png`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
#for f in `ls /global/scratch/users/drewhart/ch2/output/output/fig_hist_*png`; do rclone copy $f bdrive:ch2_outputs/analysis/; done
