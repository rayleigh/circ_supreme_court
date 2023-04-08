This repository contains the code needed to sample the dynamic circular spatial utilities model and make the WAIC comparison plots, Spearman correlation plots, and the circular evolution plots.

To run the analysis, look in the code/analysis_code/ and run run_sample_circ_id_rcpp.R after setting a result file name.

To make the plots, use the code in code/plot_code/. The results in the analysis can be used by the code in make_circ_evo_plots.R to make the circular evolution plots. For the WAIC comparison and Spearman correlation plots, a set of samples from Martin-Quinn are needed. With samples from both models, run make_waic_plots.R for the WAIC plots and make_spearman_correlation_plots.R for the Spearman correlation plots.
