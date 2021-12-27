def estimate_von_mises_mixture_params(f):
    import os
    os.system("module load r")
    os.system("Rscript --vanilla /global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/estimate_von_mises_mixture_params.R %s" % f)
    return
