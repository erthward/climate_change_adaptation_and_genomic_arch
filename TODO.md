# NEXT STEPS:

## SIMS
- debug so that each iteration gets saved (rather than only 1 iteration from each scenario,PID combo!) (but not so important because I can just run a larger number of jobs to generate the full desired number of its)

- just realized! should my null be what it currently is (i.e., no natural selection at all), or should it feature natural selection but no environmental change event??

- once everything is set, run full set of iterations, combine all output in one dir, then back up on BDrive


## ANALYSIS
- write one overarching R stats script to run all linear models on all output summary tables
- debug the skewed expectation lines in my phenotypic-space plots (only occur on Savio)

- perhaps just scatterplots instead of heatmaps for pheno shift? and then fit and plot trend line against expected, and run statistical test of area btwn trend and expected as function of genicity and linkage?

- decide if/how to run stats tests, then wrap all of that into the ch2_analys_job.sh script

- confer with Ian about final set of figs and all their details/tweaks

- once everything is set and full set of iterations are run, rerun analysis and back up full results on BDrive


## WRITING
- go back and do more reading, gather citations, and build up intro and discussion material
- begin drafting: 
  - put all results into draft as figs
  - brainstorm bulleted list of all results of interest to be mentioned
  - then freewrite rough draft of results and of discussion following those bullets


## OTHER
- consider tweaks/expansions to analysis:
  - run sims without movement (and/or with movement and dispersal restricted so that it's on contour?), to compare and show relative importance of gene flow on adaptation?
  - run sims with effect sizes drawn from distributions, to show effect of that?


# QUESTIONS FOR IAN:
- null issues:
  - null shows positive population size change for 20 but not 4 or 100 loci...
  - null have selection but no climate change (instead of no selection at all)?

- how much to toy around with running the model longer post-change?

- how much (if at all) do I need to explore other values for nuisance params?

- need to run null test for phenotypic shift, too?

- worth continuing to mess around with pop density? --> extinction happening on E edge of landsacpe in high-linkage, high-genicity scenario
  - test some metric?

- should my stats tests make linkage and genicity numeric (using recomb rates and nums of loci) or simply categorical (with pairwise t-tests to determine which categories sig diff? this allows for detection of threshold behavior, rather than assuming linear response...)?

- set of tests and outputs seem like a solid end-to-end package?

- present a summary table of my overall findings, matching my hypothesis table?
