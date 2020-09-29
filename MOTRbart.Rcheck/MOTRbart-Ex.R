pkgname <- "MOTRbart"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "MOTRbart-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('MOTRbart')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("motr_bart")
### * motr_bart

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: motr_bart
### Title: Bayesian Additive Regression Trees with Model Trees
### Aliases: motr_bart motr_bart_class

### ** Examples


# Simulate a Friedman data set
friedman_data = function(n, num_cov, sd_error){
  x = matrix(runif(n*num_cov),n,num_cov)
  y = 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5] + rnorm(n, sd=sd_error)
  return(list(y = y,
              x = x))
}
data = friedman_data(200, 10, 1)
y = data$y
x = data$x

# Run MOTR-BART
set.seed(99)
fit.motr.bart = motr_bart(x, y, ntrees = 10, nburn = 100, npost = 100)
yhat = apply(fit.motr.bart$y_hat, 2, mean)
plot(y, yhat); abline(0, 1)

# Run MOTR-BART for classification
set.seed(01)
y = ifelse(y > median(y), 1, 0)
fit.motr.bart = motr_bart_class(x, y, ntrees = 10, nburn = 100, npost = 100)
yhat = apply(fit.motr.bart$y_hat, 2, mean)
cor(y, yhat)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("motr_bart", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict_motr_bart")
### * predict_motr_bart

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict_motr_bart
### Title: Predictions for a new data set for MOTR-BART
### Aliases: predict_motr_bart predict_motr_bart_class

### ** Examples

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
y.test = ifelse(y.test > median(y.test), 1, 0)
fit.motr.bart = motr_bart_class(x, y, ntrees = 10, nburn = 100, npost = 100)
y.test.hat = predict_motr_bart_class(fit.motr.bart, x.test, 'mean')
cor(y.test, y.test.hat)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict_motr_bart", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
