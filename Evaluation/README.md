Evaluation experiments with synthetic data

The contents of this directory include the following scripts:
* `generate_synthetic_datasets.R` This file is used for generating synthetic datasets. Within this script, you can customize the seed parameters, action values, the number of classes, and the number of stages to create a generative model. To replicate the experiments conducted in the paper, please use seeds 2, 3, 4, 5, and 6.
* `results_likelihood.R` Here, the likelihood of both the training and test sets in the learned and original models is evaluated.
* `plot_results_union.R` This script displays the results and generates the convergence figure illustrating how the learned model approaches the original model.

These scripts collectively facilitate the evaluation of experiments conducted with synthetic data in the context of our research.
