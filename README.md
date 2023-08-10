# Overview

This analysis was presented at the 2023 American Ornithological Society conference, in a symposium title "Analysis challenges for breeding bird atlases".

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
2) 

