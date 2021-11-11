library(ggplot2)
library(tidyverse)
library(plyr)

# TODO:
# 1. figure out test for circular-distributed response var
# 2. decide if I need to subsample each scenario so that I have a balanced number of random chromosomes samples from the smaller-pop simulations
# 3. get running in a loop across multiple input files (for now, I can dump heads of those files into smaller example files; eventually can run on savio on full output)
# 4. figure out if/how to incorporate the enormous amount of output data into a single hierarchical model

# collect all dfs' data into one big CSV...
csvs = list.files()[grep("\\d\\.csv$", list.files())] 
dfs = list()
for (csv in csvs){
    dfs[[csv]] = read.csv(csv)
}
analysis_df = ldply(dfs, rbind)

# collect all DIR dfs' data into one big CSV...
dir_csvs = list.files()[grep("DIR_short\\.csv$", list.files())] 
#dir_csvs = list.files()[grep("DIR\\.csv$", list.files())] 
dir_dfs = list()
for (dir_csv in dir_csvs){
    dir_dfs[[dir_csv]] = read.csv(dir_csv)
}
dir_analysis_df = ldply(dir_dfs, rbind)

# collect all DIST dfs' data into one big CSV...
dist_csvs = list.files()[grep("DIST_short\\.csv$", list.files())] 
#dist_csvs = list.files()[grep("DIST\\.csv$", list.files())] 
dist_dfs = list()
for (dist_csv in dist_csvs){
    dist_dfs[[dist_csv]] = read.csv(dist_csv)
}
dist_analysis_df = ldply(dist_dfs, rbind)




# hists
# TODO: rarefy the larger-genicity data so that I've got equal random subsamples of data?
#       (does this also skew things though because it creates artificially reduced linkage
#        between separate loci included in the hist for the larger-linkage scenarios?
#        does that even matter?)
par(mfrow=c(3,3))
for (l in unique(dir_analysis_df$linkage)){
    for (g in unique(dir_analysis_df$genicity)){
        hist(dir_analysis_df[dir_analysis_df$nullness == 'non-null' & 
                             dir_analysis_df$neutrality == 'nonneut' &
                             dir_analysis_df$genicity == g &
                             dir_analysis_df$linkage == l, 'dir'],
                         breaks=500,
                         main=paste0('L=', l, '; G=', g))
    }
}

par(mfrow=c(3,3))
for (l in unique(dist_analysis_df$linkage)){
    for (g in unique(dist_analysis_df$genicity)){
        hist(dist_analysis_df[dist_analysis_df$nullness == 'non-null' & 
                             dist_analysis_df$neutrality == 'nonneut' &
                             dist_analysis_df$genicity == g &
                             dist_analysis_df$linkage == l, 'dist'],
                         breaks=500,
                         main=paste0('L=', l, '; G=', g))
    }
}




# TODO: test whether mean gene-flow vector direction is concentrated toward east in non-null sims for non-neutral loci vs neutral loci and vs non-null sims
        # TODO: use the fitted kappa vals for this??

# if so...
# for non-neutral loci, is mean direction of gene flow a fn of linkage, genicity, and their interaction?

# TODO: HOW BEST TO QUANTIZE THE LINKAGE LEVELS?
#       don't want to just assume linearity
analysis_df$linkage <- mapvalues(analysis_df$linkage, 
        from=c("independent","weak","strong"),
        to=c(0, 1, 2))

# overall dir mod
#dir_mod = lm(kappa_dir_nonneut ~ linkage * genicity,
             #data = analysis_df[analysis_df['nullness'] == 'non_null', ])
dir_mod = lm(dir ~ genicity*linkage*nullness, data=dir_analysis_df)
summary(dir_mod)
# null vs non-null mods
dir_mod_null = lm(dir ~ genicity*linkage, data=dir_analysis_df[dir_analysis_df$nullness=='null',])
summary(dir_mod_null)
dir_mod_nonnull = lm(dir ~ genicity*linkage, data=dir_analysis_df[dir_analysis_df$nullness=='non-null',])
summary(dir_mod_nonnull)


# overall dist mod
#dist_mod = lm(kappa_dist_nonneut ~ linkage * genicity,
             #data = analysis_df[analysis_df['nullness'] == 'non_null', ])
dist_mod = lm(dist ~ genicity*linkage*nullness, data=dist_analysis_df)
summary(dist_mod)
# null vs non-null mods
dist_mod_null = lm(dist ~ genicity*linkage, data=dist_analysis_df[dist_analysis_df$nullness=='null',])
summary(dist_mod_null)
dist_mod_nonnull = lm(dist ~ genicity*linkage, data=dist_analysis_df[dist_analysis_df$nullness=='non-null',])
summary(dist_mod_nonnull)



#####################################################

# vs null models...
dir_mod_null = lm(kappa_dir_nonneut ~ linkage * genicity,
             data = analysis_df[analysis_df['nullness'] == 'null', ])
summary(dir_mod_null)

# is delta_Nt and delta_fit a fn of linkage, genicity, and their interaction?
fit_mod = lm(delta_fit ~ linkage * genicity,
             data = analysis_df[analysis_df['nullness'] == 'non_null', ])
summary(fit_mod)
Nt_mod = lm(delta_Nt ~ linkage * genicity,
             data = analysis_df[analysis_df['nullness'] == 'non_null', ])
summary(Nt_mod)

# vs null models...
fit_mod_null = lm(delta_fit ~ linkage * genicity,
             data = analysis_df[analysis_df['nullness'] == 'null', ])
summary(fit_mod_null)
Nt_mod_null = lm(delta_Nt ~ linkage * genicity,
             data = analysis_df[analysis_df['nullness'] == 'null', ])
summary(Nt_mod_null)


