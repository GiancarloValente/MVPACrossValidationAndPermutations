# README	

This repository contains the accompanying code to the manuscript *Cross-validation and permutations in MVPA: validity of permutation strategies and power of cross-validation schemes* from G. Valente, A. Lage Castellanos, L. Hausfeld, F. De Martino and E. Formisano.


## Content

### Simulations

- Simulations on validity (section 3.1) and power (section 3.2.1): `SimulationsRepetitionCrossValidation`, together with the results and code to display the results.

- Decomposition of the error covariance matrix (under H0) shown in the supporting information:`ErrorCovarianceDecomposition`. The results file is rather large (1.6 GB) and not in the repository, it can be generated using the function `SH1_RunVariability_LargeNumberOfRepetitions_ErrorStats.m`.

- Comparison of unrestricted permutation and restricted permutations/randomizations with non-exchangeable data(Appendix A): `SimulationNonExchangeableData`, together with the code to display the results. When running the simulations, make sure the current directory is this folder, as some functions are redefined here and they shadow the original ones/

- Comparison of error distributions under null and alternative hypothesis (Appendix B): `ErrorDistributionComparison`. Note that the estimate of the null distribution is based on 1e6 permutations (not included in the repository)

### WU-Minn 3 Tesla dataset 

The code to prepare the data, run the analyses and display the results on the 3 Tesla WU-Minn dataset from the Human Connectome (section 3.2.2) is in the folder `ConnectomeDataset`. The data are not in the repository, to replicate the analysis the data of the 889 subjects used in this work should be downloaded from the [Human Connectome Project Database](https://db.humanconnectome.org).

## Dependencies

The code relies on the [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) implementation for the SVM and the [LibLinear](https://www.csie.ntu.edu.tw/~cjlin/liblinear/) implementation for L2-regularized logistic regression. Download and install them prior to running the current analyses, making sure their folders are added to the path.

For the analyses of the Human Connectome Project data,  a [cifti reader](https://github.com/Washington-University/cifti-matlab) should be added to the path.