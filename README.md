# Effect of income on gasoline consumption

Estimates the income elasticity of gas demand (average derivative of log gas
with respect to log income) using automatic debiased machine learning
(Riesz representer, lasso, 5-fold cross-fitting).

## Run
1. Put `gasoline_final_tf1.dta` in the same folder as the scripts.
2. In R: `source("main.R")`

## Packages
foreign, dplyr, ggplot2, quantreg, MASS, glmnet, grplasso, nnet,
randomForest, gglasso, plotrix, gridExtra, SparseM

## Output
For the full sample and each income quintile: n, AD (income elasticity), SE, 95% CI.
Full sample result: AD = 0.256, 95% CI [0.226, 0.286].
