# NEXT STEPS:

## SIMS
- debug so that each iteration gets saved (rather than only 1 iteration from each scenario,PID combo!) (but not so important because I can just run a larger number of jobs to generate the full desired number of its)

- just realized! should my null be what it currently is (i.e., no natural selection at all), or should it feature natural selection but no environmental change event??
  - right now is yes selection, no gradient, no env change --> should be yes selection, yes gradient, no change? or no to all?
  - OR SHOULD MY NULL ONLY HAVE 1 TRAIT??? whole different type of null, but in some sense this sits underneath all my hypotheses...

- once everything is set, run full set of iterations, combine all output in one dir, then back up on BDrive


## ANALYSIS
- write one overarching R stats script to run all linear models on all output summary tables
- debug the skewed expectation lines in my phenotypic-space plots (only occur on Savio)

- perhaps just scatterplots instead of heatmaps for pheno shift? and then fit and plot trend line against expected, and run statistical test of area btwn trend and expected as function of genicity and linkage?

- decide if/how to run stats tests, then wrap all of that into the ch2_analys_job.sh script

- confer with Ian about final set of figs and all their details/tweaks

- once everything is set and full set of iterations are run, rerun analysis and back up full results on BDrive


## WRITING
- then freewrite rough draft of results and of discussion following those bullets
- move gene flow directionality to last rows in the hypothesis table, not first



## OTHER
- consider tweaks/expansions to analysis:
  - run sims without movement (and/or with movement and dispersal restricted so that it's on contour?), to compare and show relative importance of gene flow on adaptation?
  - run sims with effect sizes drawn from distributions, to show effect of that?


# QUESTIONS FOR IAN:
### OVERALL

- set of tests and outputs seem like a solid end-to-end package?

### SIM Qs

  - should I try to further lengthen the 'neutral selection burn-in' period, to try to get rid of the non-zero delta_Nt and delta_fit signal in null 20-locus models (i.e., bc they're an indication that it's not hitting equilibrium prior to env change)?

##### NULL
  - null shows positive population size change for 20 but not 4 or 100 loci...
  - null have selection but no climate change (instead of no selection at all)?

### ANALYSIS Qs

  - should I write code to calculate, save, and then analyze shift in the genetic covar matrix (G) (a la McGuigan 2006 MolEcol)?
    - genetic covar matrix could be calculable, but werid because a.) usually calculated from phenotype data as a quant genetics methods casting light on allelic conditions, rather than from the genomic data themselves, and pointless because b.) I only have two traits, and I know the additive genetic variation for each (i.e., the diagonal), so it would only amount to calculating the phenotypic correlation of my two traits, which I have fixed by way of the landscape structure. What I’ve been talking about all this time is actually not a ‘shifting genetic covariance structure’, if I’m thinking about this right, but rather ‘shifting allelic covariance’, i.e., selection for new haplotypes. The sim setup feels weirder than ever now because a.) this just happens at varying rates across the landscape (i.e., env isn’t changing at the far left, then change rate increases to its maximum and far right), yet b.) at the same time gene flow can happen across the landscape, potentially swamping that signal.  ⁠⁠—> Is there a way to readjust the landscape setup so that rate of env change is constant across it?? Or should I just leave it this way and treat that as my purposeful setup. I’m back at square 1 from two years ago… :( 

  - should I consider how to look at modes of adaptation a la Retini et al. 2019, so that I can compare and contrast results to theirs?

  - how much to toy around with running the model longer post-change?

  - how much (if at all) do I need to explore other values for nuisance params?

  - need to run null test for phenotypic shift, too?

  - worth continuing to mess around with pop density? --> extinction happening on E edge of landsacpe in high-linkage, high-genicity scenario
  - test some metric?

  - should my stats tests make linkage and genicity numeric (using recomb rates and nums of loci) or simply categorical (with pairwise t-tests to determine which categories sig diff? this allows for detection of threshold behavior, rather than assuming linear response...)?

### WRITING Qs
  - include some stylized depiction/visualization of the simulation setup, with cardinal directions of gene flow labeled (N & S = ‘on contour’, E = ‘upslope’)? and include depiction of rate of env change, since it’s not homogeneous?

  - present a summary table of my overall findings, matching my hypothesis table?
