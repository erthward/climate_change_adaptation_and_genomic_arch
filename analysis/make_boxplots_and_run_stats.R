library(ggplot2)
library(tidyverse)
library(plyr)
library(nlme)
library(cowplot)
library(latex2exp)
library(ggsignif)


# TODO
# possible to add significance-testing stars above? that would be a lot of unnecessary pairwise comparisons though...


######################
# ANALYSIS FOR Nt, fit
######################

#data.dir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
data.dir = '/global/scratch/users/drewhart/ch2/output/output/'

#output.dir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
output.dir = '/global/scratch/users/drewhart/ch2/output/analysis/'

# gather all summary-output files into one
summary.csvs = list.files(data.dir)[grep("\\d\\.csv$", list.files())] 
dfs = list()
for (csv in summary.csvs){
    dfs[[csv]] = read.csv(csv)
}
summary.df = ldply(dfs, rbind)
summary.df$size = 1+0.25*floor(log10(as.numeric(as.character(summary.df$genicity))))/2
summary.df$genicity = as.factor(summary.df$genicity)

num.linkage = c()
for (linkage in summary.df$linkage){
    if (linkage == 'independent'){
        num = 0.5}
    else if (linkage == 'weak'){
        num = 0.05
    }
    else {
        num = 0.005
    }
    num.linkage = c(num.linkage, num)
}
summary.df$num.linkage = as.character(num.linkage)

plot_cols = c('#ff9a80', # lite red
              '#ed5a34', # med red
              '#c92900', # dark red
              '#80ceff', # lite blue
              '#369fe0', # med blue
              '#0368a6') # dark blue

# add color col
cols = c()
for (i in seq(nrow(summary.df))){
    nullness = summary.df[i, 'nullness']
    genicity = summary.df[i, 'genicity']
    if (genicity == 4){
        if (nullness == 'null'){
            idx = 6}
        else {idx = 3}
    }
    if (genicity == 20){
        if (nullness == 'null'){
            idx = 5}
        else {idx = 2}
    }
    if (genicity == 100){
        if (nullness == 'null'){
            idx = 4}
        else {idx = 1}
    }
    cols = c(cols, idx)
}
summary.df$color = as.factor(cols)


df.nonull = summary.df[summary.df$nullness == 'non_null',]
df.null = summary.df[summary.df$nullness == 'null',]

# plot delta_Nt and delta_fit as function of genicity and linkage
ggplot.delta_Nt = ggplot() + 
    geom_point(aes(num.linkage, delta_Nt, col=color, size=size), data=df.nonull, position='jitter') +
    geom_point(aes(num.linkage, delta_Nt, col=color, size=size), data=df.null, position='jitter') +
    scale_colour_manual(values=plot_cols) 
ggplot.delta_Nt


# boxplot instead
theme_set(theme_linedraw(base_size=20))
jpeg(paste0(output.dir, 'boxplot_delta_Nt.jpg'), width=5000, height=2500, res=300)
boxnull = ggplot(df.null) + geom_boxplot(aes(x=genicity, y=delta_Nt, fill=num.linkage)) + 
    geom_hline(yintercept=0) +
    scale_fill_manual(values = plot_cols[6:4], labels=c('strong', 'weak', 'independent')) +
    scale_y_continuous(limits = c(-750, 250)) +
    labs(y=TeX('$\\Delta$ population size'), x='number of loci per trait')

boxnonull = ggplot(df.nonull) + geom_boxplot(aes(x=genicity, y=delta_Nt, fill=num.linkage)) + 
    geom_hline(yintercept=0) +
    scale_fill_manual(values = plot_cols[3:1], labels=c('strong', 'weak', 'independent')) +
    scale_y_continuous(limits = c(-750, 250)) +
    labs(y=TeX('$\\Delta$ population size'), x='number of loci per trait')
cowplot::plot_grid(boxnull, boxnonull) 
dev.off()

# boxplot instead
jpeg(paste0(output.dir, 'boxplot_delta_fit.jpg'), width=5000, height=2500, res=300)
boxnull = ggplot(df.null) + geom_boxplot(aes(x=genicity, y=delta_fit, fill=num.linkage)) + 
    geom_hline(yintercept=0) +
    scale_fill_manual(values = plot_cols[6:4], labels=c('strong', 'weak', 'independent')) +
    scale_y_continuous(limits = c(-0.03, 0.015)) +
    labs(y=TeX('$\\Delta$ fitness'), x='number of loci per trait')
boxnonull = ggplot(df.nonull) + geom_boxplot(aes(x=genicity, y=delta_fit, fill=num.linkage)) + 
    geom_hline(yintercept=0) +
    scale_fill_manual(values = plot_cols[3:1], labels=c('strong', 'weak', 'independent')) +
    scale_y_continuous(limits = c(-0.03, 0.015)) +
    labs(y=TeX('$\\Delta$ fitness'), x='number of loci per trait')
cowplot::plot_grid(boxnull, boxnonull)
dev.off()




#######################
# ANALYSIS FOR dir DATA
#######################

# NOTE: THIS NEEDS TO HAVE BEEN RUN THROUGH THE COLUMN-SORTING PYTHON SCRIPT
#df = read.csv('./TEST_output_SORTED.csv')
#
## will run separate null and non-null tests
#df.null = df[df$nullness == 'null',]
#df.nonull = df[df$nullness == 'non-null',]
#
#
#mod.man.null = manova(cbind(mu.1, mu.2, mu.3, mu.4,
#                       kappa.1, kappa.2, kappa.3, kappa.4,
#                       alpha.1, alpha.2, alpha.3, alpha.4) ~ linkage * genicity,
#                data = df.null) 
#mod.man = manova(cbind(mu.1, mu.2, mu.3, mu.4,
#                       kappa.1, kappa.2, kappa.3, kappa.4,
#                       alpha.1, alpha.2, alpha.3, alpha.4) ~ linkage * genicity,
#                data = df.nonull) 
#
#summary.aov(mod.man.null)
#summary.aov(mod.man)



#######
# CRUFT
#######

#ggplot.delta_fit = ggplot() + 
#    geom_point(aes(genicity, delta_fit, shape=linkage), data=df.nonull, col='red') +
#    geom_point(aes(genicity, delta_fit, shape=linkage), data=df.null, col='blue')
#cowplot::plot_grid(ggplot.delta_Nt, ggplot.delta_fit)
#
#
#
## anova of delta_Nt and delta_fit as a function of genicity and linkage,
## with iteration as random fx
#mod.Nt = lme(delta_Nt ~ genicity * linkage, random = ~ 1 | .id,
#             data = summary.df[summary.df$nullness == 'non_null',])
#mod.Nt.null = lme(delta_Nt ~ genicity * linkage, random = ~ 1 | .id,
#             data = summary.df[summary.df$nullness == 'null',])
#summary(mod.Nt)
#summary(mod.Nt.null)
#
#
#
## TODO:
## 1. figure out test for circular-distributed response var
## 2. decide if I need to subsample each scenario so that I have a balanced number of random chromosomes samples from the smaller-pop simulations
## 3. get running in a loop across multiple input files (for now, I can dump heads of those files into smaller example files; eventually can run on savio on full output)
## 4. figure out if/how to incorporate the enormous amount of output data into a single hierarchical model
#
## collect all dfs' data into one big CSV...
#csvs = list.files()[grep("\\d\\.csv$", list.files())] 
#dfs = list()
#for (csv in csvs){
#    dfs[[csv]] = read.csv(csv)
#}
#analysis_df = ldply(dfs, rbind)
#
## collect all DIR dfs' data into one big CSV...
#dir_csvs = list.files()[grep("DIR_short\\.csv$", list.files())] 
##dir_csvs = list.files()[grep("DIR\\.csv$", list.files())] 
#dir_dfs = list()
#for (dir_csv in dir_csvs){
#    dir_dfs[[dir_csv]] = read.csv(dir_csv)
#}
#dir_analysis_df = ldply(dir_dfs, rbind)
#
## collect all DIST dfs' data into one big CSV...
#dist_csvs = list.files()[grep("DIST_short\\.csv$", list.files())] 
##dist_csvs = list.files()[grep("DIST\\.csv$", list.files())] 
#dist_dfs = list()
#for (dist_csv in dist_csvs){
#    dist_dfs[[dist_csv]] = read.csv(dist_csv)
#}
#dist_analysis_df = ldply(dist_dfs, rbind)
#
#
#
#
## hists
## TODO: rarefy the larger-genicity data so that I've got equal random subsamples of data?
##       (does this also skew things though because it creates artificially reduced linkage
##        between separate loci included in the hist for the larger-linkage scenarios?
##        does that even matter?)
#par(mfrow=c(3,3))
#for (l in unique(dir_analysis_df$linkage)){
#    for (g in unique(dir_analysis_df$genicity)){
#        hist(dir_analysis_df[dir_analysis_df$nullness == 'non-null' & 
#                             dir_analysis_df$neutrality == 'nonneut' &
#                             dir_analysis_df$genicity == g &
#                             dir_analysis_df$linkage == l, 'dir'],
#                         breaks=500,
#                         main=paste0('L=', l, '; G=', g))
#    }
#}
#
#par(mfrow=c(3,3))
#for (l in unique(dist_analysis_df$linkage)){
#    for (g in unique(dist_analysis_df$genicity)){
#        hist(dist_analysis_df[dist_analysis_df$nullness == 'non-null' & 
#                             dist_analysis_df$neutrality == 'nonneut' &
#                             dist_analysis_df$genicity == g &
#                             dist_analysis_df$linkage == l, 'dist'],
#                         breaks=500,
#                         main=paste0('L=', l, '; G=', g))
#    }
#}
#
#
#
#
## TODO: test whether mean gene-flow vector direction is concentrated toward east in non-null sims for non-neutral loci vs neutral loci and vs non-null sims
#        # TODO: use the fitted kappa vals for this??
#
## if so...
## for non-neutral loci, is mean direction of gene flow a fn of linkage, genicity, and their interaction?
#
## TODO: HOW BEST TO QUANTIZE THE LINKAGE LEVELS?
##       don't want to just assume linearity
#analysis_df$linkage <- mapvalues(analysis_df$linkage, 
#        from=c("independent","weak","strong"),
#        to=c(0, 1, 2))
#
## overall dir mod
##dir_mod = lm(kappa_dir_nonneut ~ linkage * genicity,
#             #data = analysis_df[analysis_df['nullness'] == 'non_null', ])
#dir_mod = lm(dir ~ genicity*linkage*nullness, data=dir_analysis_df)
#summary(dir_mod)
## null vs non-null mods
#dir_mod_null = lm(dir ~ genicity*linkage, data=dir_analysis_df[dir_analysis_df$nullness=='null',])
#summary(dir_mod_null)
#dir_mod_nonnull = lm(dir ~ genicity*linkage, data=dir_analysis_df[dir_analysis_df$nullness=='non-null',])
#summary(dir_mod_nonnull)
#
#
## overall dist mod
##dist_mod = lm(kappa_dist_nonneut ~ linkage * genicity,
#             #data = analysis_df[analysis_df['nullness'] == 'non_null', ])
#dist_mod = lm(dist ~ genicity*linkage*nullness, data=dist_analysis_df)
#summary(dist_mod)
## null vs non-null mods
#dist_mod_null = lm(dist ~ genicity*linkage, data=dist_analysis_df[dist_analysis_df$nullness=='null',])
#summary(dist_mod_null)
#dist_mod_nonnull = lm(dist ~ genicity*linkage, data=dist_analysis_df[dist_analysis_df$nullness=='non-null',])
#summary(dist_mod_nonnull)
#
#
#
######################################################
#
## vs null models...
#dir_mod_null = lm(kappa_dir_nonneut ~ linkage * genicity,
#             data = analysis_df[analysis_df['nullness'] == 'null', ])
#summary(dir_mod_null)
#
## is delta_Nt and delta_fit a fn of linkage, genicity, and their interaction?
#fit_mod = lm(delta_fit ~ linkage * genicity,
#             data = analysis_df[analysis_df['nullness'] == 'non_null', ])
#summary(fit_mod)
#Nt_mod = lm(delta_Nt ~ linkage * genicity,
#             data = analysis_df[analysis_df['nullness'] == 'non_null', ])
#summary(Nt_mod)
#
## vs null models...
#fit_mod_null = lm(delta_fit ~ linkage * genicity,
#             data = analysis_df[analysis_df['nullness'] == 'null', ])
#summary(fit_mod_null)
#Nt_mod_null = lm(delta_Nt ~ linkage * genicity,
#             data = analysis_df[analysis_df['nullness'] == 'null', ])
#summary(Nt_mod_null)
#
#
