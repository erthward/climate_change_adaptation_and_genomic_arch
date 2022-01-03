library(ggplot2)
library(tidyverse)
library(plyr)
library(nlme)
library(cowplot)
library(latex2exp)
library(ggsignif)


if (strsplit(getwd(), '/')[[1]][2] == 'home'){
    data.dir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
    output.dir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
} else {
    data.dir = '/global/scratch/users/drewhart/ch2/output/output/'
    output.dir = '/global/scratch/users/drewhart/ch2/output/analysis/'
}

###########################################################################
# delta_Nt 
  # two ANOVAs each, one null, one non-null,
  # with post-hoc all-pairwise comparisons for non-null
#--------------------------------------------------------------------------

# gather all summary-output files into one
summary.csvs = list.files(data.dir)[grep("^output_PID-\\d+\\.csv$", list.files(data.dir))] 
dfs = list()
for (csv in summary.csvs){
    dfs[[csv]] = read.csv(paste0(data.dir, csv))
}
summary.df = ldply(dfs, rbind)
summary.df$size = 1+0.25*floor(log10(as.numeric(as.character(summary.df$genicity))))/2
summary.df$genicity = as.factor(summary.df$genicity)

fact.linkage = c()
for (linkage in summary.df$linkage){
    if (linkage == 'independent'){
        num = 0.5}
    else if (linkage == 'weak'){
        num = 0.05
    }
    else {
        num = 0.005
    }
    fact.linkage = c(fact.linkage, num)
}
summary.df$fact.linkage = as.factor(as.character(fact.linkage))

df.nonull = summary.df[summary.df$nullness == 'non_null',]
df.null = summary.df[summary.df$nullness == 'null',]


cat('\n\n')
cat('DELTA_NT')
cat('\n')
cat('-------------------------------------------------------------')

# interaction plots

jpeg(paste0(output.dir, 'delta_Nt_intxn_plot.jpg'), width=5000, height=2500, res=300)
par(mfrow=c(1,2))
with(df.null, interaction.plot(genicity, fact.linkage, delta_Nt, fun = mean,
                                     main='Intxn Plot: delta_Nt: Null'))
with(df.nonull, interaction.plot(genicity, fact.linkage, delta_Nt, fun = mean,
                                     main='Intxn Plot: delta_Nt: Non-null'))
dev.off()

# fit ANOVAs
cat('\n\n')
cat('-----null---------------------------------------------------')
cat('\n\n')
mod.null <- aov(delta_Nt ~ fact.linkage + genicity + fact.linkage:genicity,
                data = df.null)
print(summary(mod.null))
# post-hoc Tukey's HSD pairwise comparison test
print("Tukey's HSD: delta_Nt: null")
print(TukeyHSD(mod.null, which = "fact.linkage"))
print(TukeyHSD(mod.null, which = "genicity"))
print(TukeyHSD(mod.null, which = "fact.linkage:genicity"))

cat('\n\n')
cat('--non-null--------------------------------------------------')
cat('\n\n')
mod.nonull <- aov(delta_Nt ~ fact.linkage + genicity + fact.linkage:genicity,
                data = df.nonull)
print(summary(mod.nonull))
print("Tukey's HSD: delta_Nt: non-null")
print(TukeyHSD(mod.nonull, which = "fact.linkage"))
print(TukeyHSD(mod.nonull, which = "genicity"))
print(TukeyHSD(mod.nonull, which = "fact.linkage:genicity"))


###########################################################################
# delta_fit
  # two ANOVAs each, one null, one non-null,
  # with post-hoc all-pairwise comparisons for non-null
#--------------------------------------------------------------------------

cat('\n\n')
cat('DELTA_FIT')
cat('\n')
cat('-------------------------------------------------------------')


# interaction plots

jpeg(paste0(output.dir, 'delta_fit_intxn_plot.jpg'), width=5000, height=2500, res=300)
par(mfrow=c(1,2))
with(df.null, interaction.plot(genicity, fact.linkage, delta_fit, fun = mean,
                                     main='Intxn Plot: delta_fit: Null'))
with(df.nonull, interaction.plot(genicity, fact.linkage, delta_fit, fun = mean,
                                     main='Intxn Plot: delta_fit: Non-null'))
dev.off()

# fit ANOVAs
cat('\n\n')
cat('-----null---------------------------------------------------')
cat('\n\n')
mod.null <- aov(delta_fit ~ fact.linkage + genicity + fact.linkage:genicity,
                data = df.null)
print(summary(mod.null))
print("Tukey's HSD: delta_fit: null")
print(TukeyHSD(mod.null, which = "fact.linkage"))
print(TukeyHSD(mod.null, which = "genicity"))
print(TukeyHSD(mod.null, which = "fact.linkage:genicity"))

cat('\n\n')
cat('--non-null--------------------------------------------------')
cat('\n\n')
mod.nonull <- aov(delta_fit ~ fact.linkage + genicity + fact.linkage:genicity,
                data = df.nonull)
print(summary(mod.nonull))
# post-hoc Tukey's HSD pairwise comparison test
print("Tukey's HSD: delta_fit: non-null")
print(TukeyHSD(mod.nonull, which = "fact.linkage"))
print(TukeyHSD(mod.nonull, which = "genicity"))
print(TukeyHSD(mod.nonull, which = "fact.linkage:genicity"))


###########################################################################
# pheno shift/shortfall
  # one ANOVA with post-hoc all-pairwise comparisons
#--------------------------------------------------------------------------

cat('\n\n')
cat('PHENO SHORTFALL')
cat('\n')
cat('-------------------------------------------------------------')


pheno.df = read.csv(paste0(output.dir, 'phenotypic_shift_undershoot.csv'))
jpeg(paste0(output.dir, 'pheno_undershoot_intxn_plot.jpg'), width=5000, height=2500, res=300)
with(pheno.df, interaction.plot(genicity, linkage, undershoot, fun = mean,
                                     main='Intxn Plot: pheno_undershoot'))
dev.off()
# make columns factors for ANOVA
pheno.df$linkage = as.factor(as.character(pheno.df$linkage))
pheno.df$genicity = as.factor(as.character(pheno.df$genicity))
# fit ANOVAs
mod <- aov(undershoot ~ linkage + genicity + linkage:genicity,
           data = pheno.df)
print(summary(mod))
# post-hoc Tukey's HSD pairwise comparison test
print("Tukey's HSD: phenotypic undershoot")
print(TukeyHSD(mod, which = "linkage"))
print(TukeyHSD(mod, which = "genicity"))
print(TukeyHSD(mod, which = "linkage:genicity"))


###########################################################################
# pop density
  # how test diff in rasters?...
#--------------------------------------------------------------------------

###########################################################################
# E- and N/S-densities
  # 4 tests: by each nullness, and for both response vars (E, N/S),
  # with bonferroni-corrected alphas
#--------------------------------------------------------------------------


