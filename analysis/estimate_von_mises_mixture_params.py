def estimate_von_mises_mixture_params(f):
    import os
    os.system("Rscript --vanilla estimate_von_mises_mixture_params.R %s" % f)
    return
