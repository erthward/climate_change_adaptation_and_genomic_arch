library(ggplot2)
library(tidyverse)
library(plyr)
library(nlme)
library(cowplot)
library(latex2exp)
library(ggsignif)



########################
# PLOT FORMATTING PARAMS
########################

legend.text.size = 17
legend.title.size = 25
axis.label.size = 25
tick.label.size = 17


######################
# ANALYSIS FOR Nt, fit
######################

if (strsplit(getwd(), '/')[[1]][2] == 'home'){
    data.dir = '/home/deth/Desktop/tmp_ch2_outputs/final_results/'
    analysis.dir = '/home/deth/Desktop/tmp_ch2_outputs/final_results/'
} else {
    datadir.file = '/global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/outputdir.txt'
    data.dir = readChar(datadir.file, file.info(datadir.file)$size-1)
    analysisdir.file = '/global/scratch/users/drewhart/ch2/climate_change_adaptation_and_genomic_arch/analysis/analysisdir.txt'
    analysis.dir = readChar(analysisdir.file, file.info(analysisdir.file)$size-1)
    #data.dir = '/global/scratch/users/drewhart/ch2/output/output/'
    #analysis.dir = '/global/scratch/users/drewhart/ch2/output/analysis/'
}

# gather all summary-output files into one
summary.csvs = list.files(data.dir)[grep("^output_PID-\\d+\\.csv$", list.files(data.dir))] 
dfs = list()
for (csv in summary.csvs){
    dfs[[csv]] = read.csv(paste0(data.dir, csv))
}
summary.df = ldply(dfs, rbind)
summary.df$size = 1+0.25*floor(log10(as.numeric(as.character(summary.df$genicity))))/2
summary.df$genicity = as.factor(summary.df$genicity)
# add row to indicate redundancy
summary.df$redundancy = ifelse(summary.df$genicity %in% c(4, 20, 100), 'lo', 'hi')

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
    if (genicity %in% c(4, 8)){
        if (nullness == 'null'){
            idx = 6}
        else {idx = 3}
    }
    if (genicity %in% c(20, 40)){
        if (nullness == 'null'){
            idx = 5}
        else {idx = 2}
    }
    if (genicity %in% c(100, 200)){
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
#ggplot.delta_Nt = ggplot() + 
#    geom_point(aes(num.linkage, delta_Nt, col=color, size=size), data=df.nonull, position='jitter') +
#    geom_point(aes(num.linkage, delta_Nt, col=color, size=size), data=df.null, position='jitter') +
#    scale_colour_manual(values=plot_cols) 
#ggplot.delta_Nt


# pop size boxplots
for (redundancy in c('lo', 'hi')){

   subdf.null = df.null[df.null$redundancy == redundancy, ]
   subdf.nonull = df.nonull[df.nonull$redundancy == redundancy, ]

  theme_set(theme_linedraw(base_size=20))
  #jpeg(paste0(analysis.dir, 'boxplot_delta_Nt_', redundancy, 'REDUND.jpg'), width=5000, height=2500, res=300)
  boxnull = ggplot(subdf.null) + geom_boxplot(aes(x=genicity, y=delta_Nt, fill=num.linkage)) + 
      geom_hline(yintercept=0) +
      scale_fill_manual(values = plot_cols[6:4],
                        labels=c('strong', 'weak', 'independent'),
                        name='linkage') +
      scale_y_continuous(limits = c(-750, 250)) +
      labs(y=TeX('$\\Delta$ population size'), x='number of loci per trait') +
      theme(legend.text=element_text(size=legend.text.size),
            legend.title=element_text(size=legend.title.size),
            axis.title=element_text(size=axis.label.size),
            axis.text=element_text(size=tick.label.size))
  
  boxnonull = ggplot(subdf.nonull) + geom_boxplot(aes(x=genicity, y=delta_Nt, fill=num.linkage)) + 
      geom_hline(yintercept=0) +
      scale_fill_manual(values = plot_cols[3:1],
                        labels=c('strong', 'weak', 'independent'),
                        name='linkage') +
      scale_y_continuous(limits = c(-750, 250)) +
      labs(y=TeX('$\\Delta$ population size'), x='number of loci per trait') +
      theme(legend.text=element_text(size=legend.text.size),
            legend.title=element_text(size=legend.title.size),
            axis.title=element_text(size=axis.label.size),
            axis.text=element_text(size=tick.label.size))
  cowplot::plot_grid(boxnull, boxnonull) 
  ggsave(paste0(analysis.dir, 'boxplot_delta_Nt_', redundancy, 'REDUND.jpg'), width=5000, height=2500, units='px', dpi=300)
  #dev.off()
  
  # mean fitness boxplots
  #jpeg(paste0(analysis.dir, 'boxplot_delta_fit_', redundancy, 'REDUND.jpg'), width=5000, height=2500, res=300)
  boxnull = ggplot(subdf.null) + geom_boxplot(aes(x=genicity, y=delta_fit, fill=num.linkage)) + 
      geom_hline(yintercept=0) +
      scale_fill_manual(values = plot_cols[6:4],
                        labels=c('strong', 'weak', 'independent'),
                        name='linkage') +
      scale_y_continuous(limits = c(-0.03, 0.015)) +
      labs(y=TeX('$\\Delta$ fitness'), x='number of loci per trait') +
      theme(legend.text=element_text(size=legend.text.size),
            legend.title=element_text(size=legend.title.size),
            axis.title=element_text(size=axis.label.size),
            axis.text=element_text(size=tick.label.size))
  boxnonull = ggplot(subdf.nonull) + geom_boxplot(aes(x=genicity, y=delta_fit, fill=num.linkage)) + 
      geom_hline(yintercept=0) +
      scale_fill_manual(values = plot_cols[3:1],
                        labels=c('strong', 'weak', 'independent'),
                        name='linkage') +
      scale_y_continuous(limits = c(-0.03, 0.015)) +
      labs(y=TeX('$\\Delta$ fitness'), x='number of loci per trait') +
      theme(legend.text=element_text(size=legend.text.size),
            legend.title=element_text(size=legend.title.size),
            axis.title=element_text(size=axis.label.size),
            axis.text=element_text(size=tick.label.size))
  cowplot::plot_grid(boxnull, boxnonull)
  ggsave(paste0(analysis.dir, 'boxplot_delta_fit_', redundancy, 'REDUND.jpg'), width=5000, height=2500, units='px', dpi=300)
  #dev.off()
}
