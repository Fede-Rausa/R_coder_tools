# R_coder_tools
some useful r functions. I have divided them in these sections:

- **regression_models**:
  includes many models for single output regression tasks. Some are invented by me (MeanSoftmaxReg and boolboost) and sometimes they works, but probably they have to be implemented in a better way. The folder ensemble_models contains three functions for bagging, boosting and bagging+boosting, that can be used for any regression model in this folder. linear_models are simple wrapper for different usages of lm function in R. gaussian_process is a really simple implementation of gaussian process, that can handle multivariate inputs and predict both mean and variance. However the implementation is smart in R, the computational cost of the model remains considerable.

- **optimization**:
  includes algorithms for free loss optimization, that I find useful in many common tasks: PSO, Simulated Annealing and Gradient Descent with numerical approximation of the gradient. All implemented from scratch.

- **preprocess**:
  includes simple functions for observing a tabular dataframe and analyze the variables inside it, plus some univariate to multivariate transformations.

- **classification_models**:
  contains some of the simplest classification models.
