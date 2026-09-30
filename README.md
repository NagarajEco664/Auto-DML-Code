# Effect of income on gasoline consumption

## project

I have replicated the work of Professor [Victor Chernozhukov](https://www.victorchernozhukov.com/) by changing the variable of interest from **price** to **income**, to understand how the method works and what it says about a different question. The original analysis estimates how gasoline demand responds to price. Here, the same automatic debiased machine learning pipeline is used to estimate how gasoline demand responds to **income**. All credit for the original method and code structure goes to the original authors. The changes in this repository adapt it to the income question.

## Question

> By what percentage does gas consumption change when household income rises by 1%?


## Changed 

| | Original | This project |
|---|---|---|
| Variable shifted to form the derivative | Price | **Income** |
| Price | Treatment | Control (price and price squared) |
| Income | Control | **Treatment** (log income, plus income squared) |
| Subgroups | Income quintiles | Income quintiles |

Every term that involves income (income squared and its interactions with the controls) is recomputed when income is shifted up and down.

## Result

The full-sample income elasticity is **0.256** (95% CI 0.226 to 0.286). A 10% rise in income goes with roughly a 2.6% rise in gas consumption. Gas is a normal good, but demand is inelastic with respect to income.

| Group | Elasticity (AD) | SE | 95% CI |
|---|---|---|---|
| Full sample | 0.256 | 0.015 | 0.226 to 0.286 |
| Quintile 1 (lowest income) | 0.232 | 0.073 | 0.087 to 0.376 |
| Quintile 2 | 0.143 | 0.071 | 0.003 to 0.283 |
| Quintile 3 | 0.174 | 0.064 | 0.049 to 0.299 |
| Quintile 4 | 0.069 | 0.062 | -0.054 to 0.191 |
| Quintile 5 (highest income) | 0.035 | 0.065 | -0.093 to 0.163 |

The response to income falls as income rises. The richest households are statistically indistinguishable from zero. The quintile intervals overlap, so the gradient is suggestive, not firmly established.

## Files

| File | Purpose |
|---|---|
| `main.R` | Runs the full analysis (full sample plus 5 income quintiles) |
| `get_data.R` | Loads the data and builds the dictionary. Income is shifted up and down to form the derivative |
| `primitives.R` | Moment functions and the finite-difference derivative with respect to log income |
| `stage1.R` | First-stage estimators of the regression and the Riesz representer (Dantzig, lasso, random forest, neural net) |
| `stage2.R` | Cross-fitting and the debiased estimate with standard error |

## Run

1. Put `gasoline_final_tf1.dta` in the same folder as the scripts.
2. In R, set the working directory to this folder and run:

```r
source("main.R")
```

The analysis takes a few minutes. It prints the sample size, elasticity (`AD`), standard error and 95% interval for each group.

## Packages

foreign, dplyr, ggplot2, quantreg, MASS, glmnet, grplasso, nnet, randomForest, gglasso, plotrix, gridExtra, SparseM

## Settings (in `main.R`)

| Setting | Options |
|---|---|
| `spec` | 1 = dictionary with income and income squared interacted with controls; 2 = adds more interactions |
| `alpha_estimator` | 0 = Dantzig, 1 = lasso |
| `gamma_estimator` | 0 = Dantzig, 1 = lasso, 2 = random forest, 3 = neural net |
| `bias` | 0 = debiased (DML), 1 = plug-in |

Avoid `gamma_estimator = 2`: a random forest is a step function, and income takes only a few bracket values, so the estimated derivative can be close to zero.

## Design notes

- **Treatment:** log income. `X.up` and `X.down` shift log income by half of `delta = sd(log income)/4` in each direction.
- **Controls:** price and price squared, age and age squared, number of drivers, household size, month, province, year, urban, young single.
- **Distance driven is deliberately excluded.** Income affects gas mostly through how much people drive, so controlling for distance blocks that channel. Including it gives an elasticity of about **0.038**, which is the effect of income on gas use *holding distance driven fixed*. That is a different and much smaller quantity than the total effect reported above.
- **Quintiles:** income quintiles are formed with `dplyr::ntile`. Income takes only a few bracket values ($20,000 to $100,000), so tied incomes are split arbitrarily and the subgroup estimates are noisy.

## Interp & limits

- The estimate is an association adjusted for the listed controls. It is causal only if those controls capture everything that moves both income and gas consumption.
- It estimates an average effect per 1% change in income, based on small shifts around each household's income. It is not a forecast for a specific household or a large income change.
- The `split.default ... not a multiple of split variable` warning is harmless. It only means the sample does not divide evenly into 5 folds.
