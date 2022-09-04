library(ggplot2)
library(ggpubr)
library(plyr)
library(tidyverse)
library(nlme)
library(cowplot)
library(latex2exp)
library(stringr)


# TODO:

  # finalize tests

  # why not 3600 rows in maladaptation?





# behavioral params
###################

factors_ordinal = F


# set up directories
####################

if (strsplit(getwd(), '/')[[1]][2] == 'home'){
    data.dir = '/media/deth/SLAB/ch2/output_100_final_for_analysis_INCOMPLETE/'
    analysis.dir = '/media/deth/SLAB/ch2/analysis/'
} else {
    datadir.file = '/global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/outputdir.txt'
    data.dir = readChar(datadir.file, file.info(datadir.file)$size-1)
    analysisdir.file = '/global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/analysisdir.txt'
    analysis.dir = readChar(analysisdir.file, file.info(analysisdir.file)$size-1)
    #data.dir = '/global/scratch/users/drewhart/ch2/output/output/'
    #analysis.dir = '/global/scratch/users/drewhart/ch2/output/analysis/'
}


# load and prep datasets
########################
# gather all demographic summary-output files into one
demog.csvs = list.files(data.dir)[grep("^output_PID-\\d+\\.csv$", list.files(data.dir))] 
dfs = list()
for (csv in demog.csvs){
    curr.df = read.csv(paste0(data.dir, csv))
    dfs[[csv]] = curr.df
}
demog.df = ldply(dfs, rbind)

# gather both maladaptation summary-output files into one
maladapt.csvs = list.files(analysis.dir)[grep("phenotypic_shift_undershoot",
                                        list.files(analysis.dir))] 
dfs = list()
for (csv in maladapt.csvs){
    # grab the redundancy value
    # (NOTE: always in a fixed location in hard-coded filenames)
    filename_info = substr(csv, 29, 100)
    if (substr(filename_info, 1, 4) == 'NULL'){
       nullness = 'null'
    } else {nullness = 'non-null'}
    if (str_count(csv, 'loREDUND') == 1){
       redundancy = 'lo'
    }
    else if (str_count(csv, 'hiREDUND') == 1){
       redundancy = 'hi'
    }
    curr.df = read.csv(paste0(analysis.dir, '/', csv))
    curr.df$nullness = nullness
    curr.df$redundancy = redundancy
    dfs[[csv]] = curr.df
}
maladapt.df = ldply(dfs, rbind)


# gather all four gene flow summary-output files into one
gf.csvs = list.files(analysis.dir)[grep("ch2_E_NS_gene_flow_densities_",
                                        list.files(analysis.dir))] 
dfs = list()
for (csv in gf.csvs){
    curr.df = read.csv(paste0(analysis.dir, '/', csv))
    # grab the redundancy value
    # (NOTE: always in a fixed location in hard-coded filenames)
    curr.df$redundancy = substr(csv, 30, 31)
    dfs[[csv]] = curr.df
}
gf.df = ldply(dfs, rbind)


# quantize independent vars' columns
quantize_ind_vars = function(df){
   # copy input df
   out_df = df

   # quantize genicity and reundancy
   stopifnot(setequal(sort(unique(out_df$genicity)), c(4, 8, 20, 40, 100, 200)))
   new_genicity = c()
   new_redundancy = c()
   for (g in out_df$genicity){
     if (g %in% c(4, 20, 100)){
        new_genicity = c(new_genicity, log(g/4, base=5)+1)
        new_redundancy = c(new_redundancy, 1)
     } else {
        new_genicity = c(new_genicity, log(g/4/2, base=5)+1)
        new_redundancy = c(new_redundancy, 2)
     }
   } 
   stopifnot(setequal(sort(unique(new_genicity)), c(1,2,3)))
   stopifnot(setequal(sort(unique(new_redundancy)), c(1,2)))
   if (factors_ordinal){
      out_df$genicity = as.numeric(new_genicity)
      out_df$redundancy = as.numeric(new_redundancy)
   } else {
      out_df$genicity = as.factor(new_genicity)
      out_df$redundancy = as.factor(new_redundancy)
   }

   # quantize linkage
   if (setequal(sort(unique(out_df$linkage)), c(0.005, 0.05, 0.5))){
      new_linkage = -log10(out_df$linkage/5)
      stopifnot(setequal(sort(new_linkage), c(1, 2, 3)))
      if (factors_ordinal){
        out_df$linkage = as.numeric(new_linkage)
      } else {
        out_df$linkage = as.factor(new_linkage)
      }
   } else if (setequal(sort(unique(out_df$linkage)), c('independent', 'strong', 'weak'))){
      new_linkage = mapvalues(out_df$linkage,
                              from=c('independent', 'weak', 'strong'),
                              to=c(1, 2, 3))
      stopifnot(setequal(sort(new_linkage), c(1, 2, 3)))
      if (factors_ordinal){
        new_linkage = as.numeric(new_linkage)
      } else {
        new_linkage = as.factor(new_linkage)
      }
   }

   # quantize nullness
   stopifnot(setequal(sort(unique(out_df$nullness)), c('non-null', 'null')) |
             setequal(sort(unique(out_df$nullness)), c('non_null', 'null')))
   new_nullness = mapvalues(out_df$nullness,
                            from=c('non-null', 'non_null', 'null'),
                            to=c(1, 1, 0))
   stopifnot(setequal(sort(unique(new_nullness)), c(0, 1)))
   if (factors_ordinal){
      out_df$nullness = as.numeric(new_nullness)
   } else {
      out_df$ nullness = as.factor(new_nullness)
   }
   return(out_df)
}
demog.df = quantize_ind_vars(demog.df)
maladapt.df = quantize_ind_vars(maladapt.df)
gf.df = quantize_ind_vars(gf.df)

# calculate derived gene flow response vars
# (as diff btwn non-null null flow densities for both onslope and upslope flow)
gf.df %>%
   group_by('genicity', 'linkage', 'redundancy') # %>%
   #mutate(E_dens_diff = 




# build and print out models
############################


# hypothesis 1.1 (upslope gene flow greater under climate change)
#   |-> expect: positive and signif coeff on nullness (coded as 0: null; 1:non-null, i.e., main)
# &
# hypothesis 1.2 (gene flow contributes least to adaptation at high polygen and low linkage)
#   |-> expect: positive and signif coeff on genicity, linkage, and their interaction term
cat('\n\n\nH1: GENE FLOW:\n------------------------------------\n\n\n')
mod.gf = lm(E_dens ~ 0 + genicity + linkage + redundancy + nullness, data=gf.df)
print(summary(mod.gf))



# hypothesis 2.1 (stronger linkage reduces adaptive capacity)
#   |-> expect: negative and signif coeff on linkage for delta_fit and delta_Nt models,
#               positive and signif coeff on linkage for maladaptation model
# &
# hypothesis 2.2 (higher polygenicity reduces adaptive capacity)
#   |-> expect: negative and signif coeff on genicity for delta_fit and delta_Nt models,
#               positive and signif coeff on genicity for maladaptation model
# &
# hypothesis 3 (greater redundancy increases adaptive capacity)
#   |-> expect: positive and signif coeff on redundancy for delta_fit and delta_Nt models,
#               ngeative and signif coeff on redundancy for maladaptation model
cat('\n\n\nH2 & H3: ADAPTIVE CAPACITY:\n------------------------------------\n\n\n')
mod.delta_fit = lm(delta_fit ~ 0 + genicity + linkage + redundancy + nullness, data=demog.df)
print(summary(mod.delta_fit))
cat('\n\n')
mod.delta_Nt = lm(delta_Nt ~ 0 + genicity + linkage + redundancy + nullness, data=demog.df)
summary(mod.delta_Nt)
mod.undershoot = lm(undershoot ~ 0 + genicity + linkage + redundancy + nullness, data=maladapt.df)
print(summary(mod.undershoot))



