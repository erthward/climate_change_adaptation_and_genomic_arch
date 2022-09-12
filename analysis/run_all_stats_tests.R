library(ggplot2)
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

# treat the independent variables as ordinal factors,
# instead of continuous numerical vars?
ind_vars_as_ordinal_factors = F


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
}


# load and prep datasets
########################
# gather all demographic and gene flow summary-output files into one
demog.gf.csvs = list.files(data.dir)[grep("^output_PID-\\d+\\.csv$", list.files(data.dir))] 
dfs = list()
for (csv in demog.gf.csvs){
    curr.df = read.csv(paste0(data.dir, csv))
    dfs[[csv]] = curr.df
}
demog.gf.df = ldply(dfs, rbind)

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
        new_genicity = c(new_genicity, log(g/4, base=5))
        new_redundancy = c(new_redundancy, 0)
     } else {
        new_genicity = c(new_genicity, log(g/4/2, base=5))
        new_redundancy = c(new_redundancy, 1)
     }
   } 
   stopifnot(setequal(sort(unique(new_genicity)), c(0,1,2)))
   stopifnot(setequal(sort(unique(new_redundancy)), c(0,1)))
   if (ind_vars_as_ordinal_factors){
      out_df$genicity = as.factor(new_genicity)
      out_df$redundancy = as.factor(new_redundancy)
   } else {
      out_df$genicity = as.numeric(new_genicity)
      out_df$redundancy = as.numeric(new_redundancy)
   }

   # quantize linkage
   if (setequal(sort(unique(out_df$linkage)), c(0.005, 0.05, 0.5))){
      new_linkage = -log10(out_df$linkage/0.5)
      stopifnot(setequal(sort(new_linkage), c(0, 1, 2)))
      if (ind_vars_as_ordinal_factors){
        out_df$linkage = as.factor(new_linkage)
      } else {
        out_df$linkage = as.numeric(new_linkage)
      }
   } else if (setequal(sort(unique(out_df$linkage)), c('independent', 'strong', 'weak'))){
      new_linkage = mapvalues(out_df$linkage,
                              from=c('independent', 'weak', 'strong'),
                              to=c(0, 1, 2))
      stopifnot(setequal(sort(new_linkage), c(0, 1, 2)))
      if (ind_vars_as_ordinal_factors){
        new_linkage = as.factor(new_linkage)
      } else {
        new_linkage = as.numeric(new_linkage)
      }
      out_df$linkage = new_linkage
   }

   # quantize nullness
   stopifnot(setequal(sort(unique(out_df$nullness)), c('non-null', 'null')) |
             setequal(sort(unique(out_df$nullness)), c('non_null', 'null')))
   new_nullness = mapvalues(out_df$nullness,
                            from=c('non-null', 'non_null', 'null'),
                            to=c(1, 1, 0))
   stopifnot(setequal(sort(unique(new_nullness)), c(0, 1)))
   # NOTE: always treat nullness as categorical
   out_df$ nullness = as.factor(new_nullness)
   return(out_df)
}
demog.gf.df = quantize_ind_vars(demog.gf.df)
maladapt.df = quantize_ind_vars(maladapt.df)

cat('\n\n\n')
cat('SAMPLE SIZES:\n-------------\n')
cat('demog.gf.df has ', nrow(demog.gf.df), ' rows\n\n')
cat('maladapt.df has ', nrow(maladapt.df), ' rows')
cat('\n\n\n')

# calculate derived gene flow response var
# (as diff btwn non-null null flow density upslope)
# NOTE: coded as 'Eness' (i.e., 'eastness')
null.gf.Eness.df = demog.gf.df %>%
   group_by('.id', 'genicity', 'linkage', 'redundancy') %>%
   subset(subset=nullness==0)
gf.Eness.df = demog.gf.df %>%
   group_by('.id', 'genicity', 'linkage', 'redundancy') %>%
   subset(subset=nullness==1)
gf.Eness.df['delta_flow'] = gf.Eness.df['Eness'] - null.gf.Eness.df['Eness']


# print structure info for each data.frame, to check that factor/numeric variables are correct
print(str(demog.gf.df))
print(str(maladapt.df))
print(str(gf.Eness.df))

# build and print out models
############################


# hypothesis 1.1 (upslope gene flow greater under climate change)
#   |-> expect: positive and signif coeff on nullness (coded as 1: null; 2:non-null, i.e., main)
# &
# hypothesis 1.2 (gene flow contributes least to adaptation at high polygen and low linkage)
#   |-> expect: positive and signif coeff on genicity, linkage, and their interaction term
cat('\n\n\nH1: GENE FLOW:\n------------------------------------\n\n\n')
mod.gf = lm(delta_flow ~ genicity + linkage + redundancy, data=gf.Eness.df)
print(summary(mod.gf))

# predict values for each scenario
pred_vals = data.frame(expand.grid(unique(demog.gf.df$genicity),
                                   unique(demog.gf.df$linkage),
                                   unique(demog.gf.df$redundancy)))
colnames(pred_vals) = c('genicity', 'linkage', 'redundancy')
pred_vals$PREDICTED_delta_flow = predict.lm(mod.gf, pred_vals)
cat('\n\tpredict values of change in mean upslope gene flow:\n')
print(pred_vals)
cat('\n\n')


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
mod.delta_fit = lm(delta_fit ~ genicity + linkage + redundancy + nullness, data=demog.gf.df)
print(summary(mod.delta_fit))
cat('\n\n')
mod.delta_Nt = lm(delta_Nt ~ genicity + linkage + redundancy + nullness, data=demog.gf.df)
summary(mod.delta_Nt)
mod.maladapt = lm(undershoot ~ genicity + linkage + redundancy + nullness, data=maladapt.df)
print(summary(mod.maladapt))
# predict values for each scenario
pred_vals = data.frame(expand.grid(unique(demog.gf.df$nullness),
                                   unique(demog.gf.df$genicity),
                                   unique(demog.gf.df$linkage),
                                   unique(demog.gf.df$redundancy)))
colnames(pred_vals) = c('nullness', 'genicity', 'linkage', 'redundancy')
pred_vals$PREDICTED_delta_fit = predict.lm(mod.delta_fit, pred_vals)
pred_vals$PREDICTED_delta_Nt = predict.lm(mod.delta_Nt, pred_vals)
pred_vals$PREDICTED_maladapt = predict.lm(mod.maladapt, pred_vals)
cat('\n\tpredicted values of change in fitness, change in Nt, and maladaptation:\n')
print(pred_vals)
cat('\n\n')


