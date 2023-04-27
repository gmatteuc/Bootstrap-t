# Bootstrap-t

Matlab implementation of Efron's Bootstrap-t procedure to estimate confidence intervals. Follows the procedure described in:
"Tibshirani, R. J., & Efron, B. (1993). An introduction to the bootstrap. Monographs on statistics and applied probability, 57(1)"

Provided any function (taking as input N M-dimensional matrices) performs Bootstrap-t by resampling on a specified dimension of each input maxtrix and returns confidence intervals for all the outputs of the function (any K L-dimensional matrices). The code is designed to be general and easily adapt to any function data for which is hard to estimate confidence intervals otherwhise. The code can run the outer bootstrap resampling in a parallelized way to speed up computation time.

This repository includes:

test_Bootstrap_t.m ---> a script to test and show an example of usage of the provided functions
get_mean_difference_Bootstrap_t.m ---> wrapper function encapsulating a mean difference (function replaceable by some user provided function) for the test script
get_Bootstrap_t_ci_serial ---> core function implementing the Bootstrap-t algorithm (non-parallel implementation)
get_Bootstrap_t_ci_parallel ---> core function implementing the Bootstrap-t algorithm (parallel implementation)
cell2matlastdim.m ---> utility function for cell array content concatenation

Code written and tested in MATLAB R2019b.
Author: Giulio Matteucci
Date: 21/04/2023
