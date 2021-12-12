# NEXT STEPS:

## SIMS
- debug KDE plot
- debug so that each iteration gets saved (rather than only 1 iteration from each scenario,PID combo!)
- improve setup:
  - need to run a 'selection burn-in' phase, until populations and average fitnesses are stable --> otherwise, current pre-climate-change phase is not serving as a proper comparison because the population still is not in equilibrium
  - do something to constrain and/or measure the lateral spread of a pop (b/c see 100-gene sims, where pop is intially in a lateral band, then spreads east... that doesn't occur anywhere else... interesting result, or confounding of results?)

## ANALYSIS
1.)
- decide if/how to run stats tests!
- wrap all of that into the ch2_analys_job.sh script

2.)
- write a short and sweet script to loop through output files and visualize population size changes and fitness changes, with uncertainty bounds
- wrap that into ch2_analysis_job.sh


## OTHER
- need to standardize y axes on histograms, pop-size plots, fitness plots?
- figure out why shift and null don't seem to be plotting together correctly
  on map and pop-size plots
