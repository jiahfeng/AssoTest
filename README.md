
**Table of Contents**

- [Package](#Package)
- [Installation](#Installation)
- [Usage](#Usage)
  * [scoretest](#scoretest)
  * [surv_score](#surv_score)
- [Help](#Help)
- [Contact](#Contact)
- [Reference](#Reference)


# Package
<ins>**AssoTest**</ins> is a package that performs score test to investigate the association between an outcome variable and an incomplete covariate, using the method proposed by [Wong et al.](https://doi.org/10.5705/ss.202021.0253) (2023). 
**AssoTest** relies on the R-packages `glmnet`, `stats`, `Rcpp` and `RcppEigen`, which are hosted on CRAN.


# Installation 
**AssoTest** can be installed from Github directly:
```
install.packages("devtools")
library(devtools)
install_github("jiahfeng/AssoTest")
library(AssoTest)
```

# Usage
The package contains 2 main functions:
Function  | Description
------------- | -------------
scoretest  | Performs score test with missing data for continuous and binary outcomes
surv_score  |  Performs (supremum) score test with missing data for time-to-event outcomes

## scoretest

```
scoretest(x, y, s, w = NULL, r, family = c("gaussian", "binomial"))

```
This function performs a score test to investigate the association between a phenotype and an incomplete covariate, where the incomplete covariate may be associated with potentially high-dimensional auxiliary variables.

The details of the model can be found in [Wong et al.](https://doi.org/10.5705/ss.202021.0253) (2023).


## surv_score

```
scoretest_supre(x, y, s, w, r, delta, Z)

```
This function performs a score test to investigate the association between a right-censored survival outcome and an incomplete covariate, where the incomplete covariate may be associated with potentially high-dimensional auxiliary variables. We specify the link between the survival outcome and the covariates using a transformation model, which includes the proportional hazards model and the proportional odds model as special cases.
 
### Example

```
load("example_data.RData")
scoretest_supre(x = data_list$x, y = data_list$y, s = data_list$s, w = data_list$w, r = data_list$r, delta = data_list$delta, Z= data_list$Z)
```

# Help 

Details about the package can be found in the user manual:
```
?scoretest
```

# Contact 
Jiahui Feng <<jhfeng731@gmail.com>>

# Reference 
Wong, K. Y., & Feng, J. (2023). Score tests with incomplete covariates and high-dimensional auxiliary variables. _Statistica Sinica_. **33**: 1483-1505.

