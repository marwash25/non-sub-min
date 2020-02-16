# non-sub-min

Matlab code for reproducing the experimental results of the paper [Optimal approximation for unconstrained non-submodular minimization](https://arxiv.org/abs/1905.12145), by Marwa El Halabi and Stefanie Jegelka.

## Usage:

### Noisy submodular minimization experiment

To reproduce the figures for this experiment run the script `noisy_sfm_plot.m` (with `dataset = 1` for "Genrmf-long", and `dataset = 3` for "Two-moons").

To reproduce the results call the function:

`noisy_sfm_test(run,dataset,ind_algo,m,sigma,maxiter_in)`

with the following parameters (see paper for details):
* `run` index of the run
*	`dataset` index of the dataset (1 - "Genrmf-long", 3- "Two-moons")
*	`ind_algo` index of the algorithm (1- MNP, 2- CG-LS, 3- CG, 4- PGM, 5- PGM-polyak, 6- ACCPM, 7- ACCPM-Kelley)
*	`m` number of samples
*	`sigma` standard deviation of noise
*	`maxiter_in` maximum number of iterations 

**Credits**: The code for this experiment is from [B16], with some modifications.

### Structured sparse learning experiment

To reproduce the figures for this experiment run the script `Range_plot.m` 

To reproduce the results call the function:

`Range_test(run,d,k,ind_n,sig_noise,corr_level,maxiter)`

with the following parameters (see paper for details):
* `run` index of the run
* `d` dimension
* `k` sparsity
* `ind_n` index specifying number of measurements (from `n_vec = ceil(linspace(0.25,2,10)*d)`)
* `sig_noise` standard deviation of noise
* `corr_level` correlation level of measurements (set to 0)
* `maxiter` maximum number of iterations for PGM

**Credits**: The implementation of PGM is a slightly modified version of the code from [B16]. 
The code to compute rank-1 updates of the pseudo-inverse of a matrix is from [K19].

**Dependencies**: CVX [GB14] 

## References:
 
* [B16] F. Bach. Submodular package, version 2.0. https://www.di.ens.fr/~fbach/submodular/, retrieved Sept. 2016.
* [K19] Navvab Kashiri. Moore-Penrose Pseudo-Inverse Rank-1 Update, MATLAB Central File Exchange. https://www.mathworks.com/matlabcentral/fileexchange/61115-moore-penrose-pseudo-inverse-rank-1-update, retrieved April 2019.
* [GB14] M. Grant and S. Boyd. CVX: Matlab software for disciplined convex programming, version 2.1.
http://cvxr.com/cvx/, Mar. 2014.


