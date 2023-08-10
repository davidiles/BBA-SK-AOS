# Overview
_Written by David Iles on 10 August 2023_

This analysis was presented at the 2023 American Ornithological Society conference, in a symposium title "Analysis challenges for breeding bird atlases".  This is mostly a personal record of the analyses I am developing for Canada's breeding bird atlases, and there is no guarantee these scripts are "final". 

# Objectives

Use 5-fold spatially blocked cross-validation to evaluate the costs/benefits of:

1) Integrated analysis of point counts and presence/absence checklists (stationary and linear transects)
2) Increased point count coverage across the province (using half of squares for training, versus using all squares for training)
3) Filling gaps between point count coverage with checklists

Cross-validation accuracy is measured using:

1) Area under the curve (ability to accurately predict species presences/absences)
2) Correlation between predicted and observed counts (ability to predict spatial patterns)
3) Mean squared error (ability to predict counts)
4) Log-pointwise predictive density (ability to predict counts)

# Workflow

1) `1_select_xval_squares.R` divides the landscape into a spatially balanced, random selection of 50 km x 50 km blocks.  It then separates them into 5 cross-validation folds.
2) `2_Prep_Bird_Data.R` takes survey data (point counts and checklists) provided by Birds Canada (through NatureCounts) and cleans/formats it for analysis
3) `3_prepare_AOS_data.R` sets the criteria for inclusion in analysis (dates of surveys, times of day, etc).  This could be merged with `2_Prep_Bird_Data.R` in future projects.
4) `4_AOS_xval_PC_CL_subset.R` conducts the integrated analyses and saves cross-validation accuracy metrics
5) `5_AOS_plots_subset.R` summarizes and plots results of cross-validation analysis.
6) `6_AOS_surface_comparison_subset.R` creates a full prediction surface across SK from each model, using all available data.  This would be the "full analysis", but some pieces of code need to be updated.

# Other notes

- the "other_analyses" folder contains scripts that conduct the crossvalidation analyses using other data selection criteria.  These scripts are a work in progress.
- the "actual" provincial breeding bird atlas analysis will be contained in a different github repository.



