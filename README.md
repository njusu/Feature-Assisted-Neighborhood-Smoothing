## Feature-Assisted-Neighborhood-Smoothing
#### FANS.R
  Input --
  A: adjacency matrix;
  X: feature matrix;
  lambda: tuning parameter;
  screen: controls if feature screening is performed.
  Output:
  P_hat: estimated graphon values or probability matrix.
#### CV_for_FANS.R
  Takes a list of lambda values (LAMBDA) as input.
  Outputs the CV error for each lambda.
  
#### Other related algorithms are in "Matlab Code" folder

## Plots in the paper
#### visualize_graphon.R
  Figure 1, 2 (top), 3, 10, 11
#### MSE_vs_lambda.R
  Figure 2 (bottom), 4
#### repeatCompare.m (in "Matlab Code" folder)
  Figure 5
#### realdata_results.R and "Matlab Code":
  Figure 6, 7, 8, 9
