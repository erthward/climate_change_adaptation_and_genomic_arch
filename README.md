# Code from: Genomic architecture controls multivariate adaptation to climate change
### (chapter 2 of Drew E. Terasaki Hart's PhD dissertation)

Code: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10531208.svg)](https://doi.org/10.5281/zenodo.10531208)

Data: 

-----------


#### Workflow used:

1. set `redundancy = 'hi'` in ./sim/batch_script.py
2. run 5 jobs of ./sim/ch2_sim_job.sh (each producing 20 sets of results) 
3. run ./sim/count_complete_datasets.py with args `f t` to check if 100 complete high-redundancy datasets available
4. if so, run with args `t t` to ensure assert statemnts pass
5. if so, run with args `t t /PATH/TO/output_100_final/` to copy all 100 datsets to that directory
6. set `redundancy = 'lo'` in ./sim/batch_script.py, then repeat steps 2-5 (but with second arg to ./sim/count_complete_datasets.py always `f` instead of `t`, to indicate low redundancy)
7. copy the full set of 100 datsets to Google Drive for backup
8. run all analysis jobs (./analysis/ch2_Nt_mean_fit_analysis_job.sh, ./analysis/ch2_pheno_shift_analysis_job.sh, ./analysis/ch2_gene_flow_analysis_job.sh, and ./analysis/ch2_stats_pop_dens_analysis_job.sh), to generate results and push them to Google Drive
9. use results to compose final figures, then move those figures to ./pub
10. compile paper and submit


