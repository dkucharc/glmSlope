# glmSLOPE - SLOPE for generalized linear models

[![Build Status](https://travis-ci.com/dkucharc/glmSLOPE.svg?branch=master)](https://travis-ci.com/dkucharc/glmSLOPE)

### Instalation guide
Currently the package is under heavy development. Therefore, it is recommended to install the package directly from the Github repository. This can achieved by using an `install_github` function provided by the `devtools` package. Namely,
```R
library(devtools)
install_github("dkucharc/glmSLOPE")
```

### Examples
- **Linear model**
```R
library(glmSLOPE)
# Test data
X <- c(0.53766714, 1.833885, -2.2588469, 0.86217332, 0.31876524, -1.3076883, -0.43359202, 0.34262447, 3.5783969, 2.769437, -1.3498869, 3.0349235, 0.72540422, -0.063054873, 0.7147429, -0.20496606, -0.12414435, 1.4896976, 1.4090345, 1.4171924, 0.67149713, -1.2074869, 0.71723865, 1.6302353, 0.48889377, 1.034693, 0.72688513, -0.30344092, 0.29387147, -0.7872828, 0.88839563, -1.1470701, -1.0688705, -0.80949869, -2.9442842, 1.4383803, 0.32519054, -0.75492832, 1.3702985, -1.7115164, -0.10224245, -0.24144704, 0.31920674, 0.3128586, -0.86487992, -0.030051296, -0.16487902, 0.62770729, 1.0932657, 1.109273)
dim(X) <- c(5, 10)
y <- c(1.0734014, -5.3021346, 1.096639, -0.39124089, -0.92884291)
dim(y) <- c(5,1)
lambda <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
# Estimate parameters
solve_slope(X, y, lambda, model = 'linear')
```

- **Logistic model**
```R
library(glmSLOPE)
# Test data
X <- c(0.53766714, 1.833885, -2.2588469, 0.86217332, 0.31876524, -1.3076883, -0.43359202, 0.34262447, 3.5783969, 2.769437, -1.3498869, 3.0349235, 0.72540422, -0.063054873, 0.7147429, -0.20496606, -0.12414435, 1.4896976, 1.4090345, 1.4171924, 0.67149713, -1.2074869, 0.71723865, 1.6302353, 0.48889377, 1.034693, 0.72688513, -0.30344092, 0.29387147, -0.7872828, 0.88839563, -1.1470701, -1.0688705, -0.80949869, -2.9442842, 1.4383803, 0.32519054, -0.75492832, 1.3702985, -1.7115164, -0.10224245, -0.24144704, 0.31920674, 0.3128586, -0.86487992, -0.030051296, -0.16487902, 0.62770729, 1.0932657, 1.109273)
dim(X) <- c(5, 10)
y <- c(1, -1, 1, -1, -1)
dim(y) <- c(5,1)
lambda <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
# Estimate parameters
solve_slope(X, y, lambda, model = 'logistic')
```

### References:
1. [SLOPE—Adaptive variable selection via convex optimization](https://projecteuclid.org/download/pdfview_1/euclid.aoas/1446488733)

2. [Sparse Portfolio Selection via the sorted L1 - Norm](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3048879&download=yes)

3. [Proximal Algorithms](https://web.stanford.edu/~boyd/papers/pdf/prox_algs.pdf)
