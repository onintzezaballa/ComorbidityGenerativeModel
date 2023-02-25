# A Probabilistic Generative Model to Discover the Treatments of Coexisting Diseases with Missing Data
[![License-MIT](https://img.shields.io/badge/License-MIT-red)](/LICENSE) [![Made-with-R](https://img.shields.io/badge/Made%20with-R-blue)](/GenerativeModel) [![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](#contact)
 
This repository contains the code developed for the paper "A Probabilistic Generative Model to Discover the Treatments of Coexisting Diseases with Missing Data". This paper presents a probabilistic generative model of coexisting treatments that are described in terms of sequences of medical activities of variable length. We consider the specific problem of learning the joint progression of comorbidities where most of the diagnoses are missing. The generative model has a three-fold objective: (i) identify and segment the medical history of patients into treatments associated with comorbidities; (ii) learn the model associated with each identified disease treatment; and (iii) discover subtypes of patients with similar coevolution of comorbidities. To this end, the model considers a latent structure for the sequences, where patients are modeled by a latent class defined by the evolution of their comorbidities, and each observed medical event of their clinical history is asso- ciated with a latent disease. The learning process is performed using an Expectation-Maximization algorithm that considers the exponential number of configurations of the latent variables and is efficiently solved with dynamic pro- gramming. The evaluation of the method is carried out both on synthetic and real world data: the experiments on synthetic data show that the learning procedure allows the generative model underlying the data to be recovered; the experiments on real medical data show accurate results in the segmentation of sequences into different treatments, subtyping of patients and diagnosis imputation.

It is a free R code that is under [MIT License](/LICENSE).

## Implementation of the method

[GenerativeModel](/GenerativeModel) folder contains R scripts required to execute the method:

* `main.R` is the main file. In this file we can modify the values of the parameters, such as the values of actions, the number of classes and the number of stages, to create a generative model. In addition, such script randonmly samples sequences of actions and learns a new model with the EM algorithm. Afterwards, the error between the learned model and the original one is computed to show that the proposed method recovers the real model unerlying the data.
* `generateModel.R` generates the original model and samples sequences of actions from such probabilistic model.
* `initialization.R` initializes the parameters of the model
* `EMalgorithm.R` efficiently performs the Expectation-Maximization algorithm, where the parameters are updated in each iteration with a dynamic programming based method.
* `smoothing.R` contains functions to carry out the smoothing and normalization of the parameters of a model.
* `Loglikelihod.R` computes the log-likelihood of a model
* `MSE.R` computes the Mean Squared Error of a trained model by comparing with the original model which the sequence were sampled from.


## Data

We use synthetic datasets generated with the proposed model.


## Evaluation

We display in this repository an evaluation of the model to demonstrate that we are able to recover the original generative model, similar to the one performed on the paper. [Evaluation](/Evaluation) folder includes more details of the dataset generation, commands to execute the method, and results. 



## Contact
Onintze Zaballa Larumbe

onintzezaballa@gmail.com

[![ForTheBadge built-with-science](http://ForTheBadge.com/images/badges/built-with-science.svg)](https://github.com/onintzezaballa)

