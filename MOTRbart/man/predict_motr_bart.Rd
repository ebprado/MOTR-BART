\name{predict_motr_bart}
\alias{predict_motr_bart}
\alias{predict_motr_bart_class}
\title{Predictions for a new data set for MOTR-BART}

\description{
This function generates predictions for a new data set from a MOTR-BART object.}

\usage{
predict_motr_bart(object,
                  newdata,
                  type = c("all", "median", "mean"))

predict_motr_bart_class(object,
                        newdata,
                        type = c("all", "median", "mean"))

}
\arguments{
\item{object}{The fit object returned by the function \code{motr_bart}.}

\item{newdata}{The design matrix that will be used in the prediction.}

\item{type}{Should the predicted values be the \code{"mean"} or \code{"median"} of the posterior draws? It is also possible to set \code{"all"}, where no function is applied and the posterior draws are returned for each MCMC iteration.}
}
\value{
If \code{type = "mean"} or \code{type = "median"}, then a vector with predicted values is returned. If \code{type = "all"}, an n times npost is returned, where n is the number of row of \code{newdata}.

Note that when the response variable is binary, \code{y_hat} contains the draws of the probability P(Y = 1 | x).

}

\examples{
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
}
