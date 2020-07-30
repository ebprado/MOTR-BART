# BART with Model Trees (MOTR-BART)

This repository contains R scripts and data sets that can be used to reproduce the results presented in "Estevão B Prado, Rafael A Moral, Andrew C Parnell (2020). Bayesian Additive Regression Trees with Model Trees. ArXiv preprint ArXiv: 2006.07493".

In addition, it provides an implementation of MOTR-BART in a format of R package named ```MOTRbart```.

## Installation
``` r
library(devtools)
install_github("ebprado/MOTR-BART/tree/master/MOTRbart")
```
## Example
``` r
library(MOTRbart)
# Simulate a Friedman data set
friedman_data = function(n, num_cov, sd_error){
  x = matrix(runif(n*num_cov),n,num_cov)
  y = 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5] + rnorm(n, sd=sd_error)
  return(list(y = y,
              x = x))
}
# Training data
data = friedman_data(200, 10, 1)
y = data$y
x = data$x

# Test data
data_test = friedman_data(100, 10, 1)
y.test = data_test$y
x.test = data_test$x

# Run MOTR-BART
set.seed(99)
fit.motr.bart = motr_bart(x, y, ntrees = 10, nburn = 100, npost = 100)
y.test.hat = predict_motr_bart(fit.motr.bart, x.test, 'mean')
plot(y.test, y.test.hat); abline(0, 1)

# Run MOTR-BART for classification
set.seed(01)
y = ifelse(y > median(y), 1, 0)
y.test = ifelse(y.test > median(y), 1, 0)
fit.motr.bart = motr_bart_class(x, y, ntrees = 10, nburn = 100, npost = 100)
y.test.hat = predict_motr_bart_class(fit.motr.bart, x.test, 'mean')
cor(y.test, y.test.hat)
```
