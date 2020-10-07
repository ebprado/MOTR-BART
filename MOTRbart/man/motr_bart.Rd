\name{motr_bart}
\alias{motr_bart}
\alias{motr_bart_class}
\title{Bayesian Additive Regression Trees with Model Trees}
\description{
MOTR-BART is a Bayesian non-parametric model where the predictions are based on piece-wise linear functions rather than terminal node constants.

\itemize{
  \item For numeric outcomes \eqn{y = f(x) + \epsilon}{y = f(x) + \epsilon},
where \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)}.
  \item For binary outcomes \eqn{y}, \eqn{P(Y = 1 \mid x) = \Phi(f(x))}{P(Y = 1 | x) = \Phi(f(x))}, where \eqn{\Phi}
denotes the cdf of the standard Normal distribution.
}
}
\usage{
motr_bart(x,
          y,
          sparse = TRUE,
          vars_inter_slope = TRUE,
          ntrees = 10,
          node_min_size = 5,
          alpha = 0.95,
          beta = 2,
          nu = 3,
          lambda = 0.1,
          sigma2 = 1,
          nburn = 1000,
          npost = 1000,
          nthin = 1,
          ancestors = FALSE)

motr_bart_class(x,
                y,
                sparse = TRUE,
                vars_inter_slope = TRUE,
                ntrees = 10,
                node_min_size = 5,
                alpha = 0.95,
                beta = 2,
                nu = 3,
                lambda = 0.1,
                sigma2 = 1,
                nburn = 1000,
                npost = 1000,
                nthin = 1,
                ancestors = FALSE)
}

\arguments{
   \item{x}{
     The design matrix containing the covariates that will be used in the regression.
   }

   \item{y}{
      The response variable (continuous or binary).
   }

   \item{sparse}{
   If \code{TRUE}, then the variable selection is performed based on the Dirichlet prior as in Linero (2016). If FALSE, then the variable selection is based on a uniform prior as in the standard BART. Default is \code{TRUE}.
   }
   \item{vars_inter_slope}{
   If \code{TRUE}, then the intercepts and the slopes are updated and penalised differently. Default is \code{TRUE}.
   }
   \item{ntrees}{
    Number of trees. Default is 10.
   }

   \item{node_min_size}{
   	 Size of the smallest node. Default is 5.
   }

   \item{alpha}{
     Parameter of the prior for the tree. It should be a value within (0,1). Default is 0.95.
   }

   \item{beta}{
     Parameter of the prior for the tree. It should be a value greater or equal to 0. Default is 2.
   }

   \item{nu}{
     Shape of the prior for the variance.
   }

   \item{lambda}{
     Scale of the prior for the variance.
   }

   \item{sigma2}{
     An inital estimate of the variance.
   }

   \item{nburn}{
      Number of MCMC iterations for the burn-in period.
  }
  \item{npost}{
     Number of MCMC iterations that should be performed after the burn-in period.
   }
  \item{nthin}{
     Amount of thinning that should be considered. If it is greater than 1, then only the draws every "nthin" iterations will be returned.
   }
   \item{ancestors}{
Should the ancestors be included in the linear predictors. Default is \code{FALSE}. If \code{ancestors = TRUE}, then only those covariates that are ancestors of the terminal nodes are used in the linear predictors. If \code{ancestors = FALSE}, then only those covariates that are used in the splitting rules of the corresponding tree are considered in the linear predictors.
}
}

\value{
The following objects are returned by the \code{motr_bart} and \code{motr_bart_class} functions. When the response variable is binary, \code{y_hat} contains the draws of the probability P(Y = 1 | x).

\item{\code{y_hat}}{
An (npost x n) matrix with the predicted values for each MCMC iteration.}

\item{\code{center_x}}{
Mean of the continuous covariates in the design matrix. For binary covariates, the mean is 0.}

\item{\code{scale_x}}{
Standard deviation of the continuous covariates in the design matrix. For binary covariates, the standard deviation is 1.}

\item{\code{npost}}{
Number of MCMC iterations after the burn-period.}

\item{\code{nburn}}{
Number of MCMC iteration as burn-in}

\item{\code{nthin}}{
Amount of thinning.}

\item{\code{ntrees}}{
Number of trees.}

\item{\code{y_mean}}{
Mean of the response variable.}

\item{\code{y_sd}}{
Standard deviation of the response variable.}

\item{var_count_store}{
A matrix with \code{npost} rows and p columns. It returns the total count of the number of times that the covariates are used over all trees decision rule.
}
\item{s}{
A matrix with \code{npost} rows and p columns. It returns the posterior distribution of s; see Linero (2016).
}

\item{vars_betas}{
A matrix with the posterior variances of the intercepts and slopes. The first column refers to the variance of the intercepts and the second one to the variance of the slopes.
}
}

\references{
Prado, E. B., Moral, R. A., & Parnell, A. C. (2020). Bayesian Additive Regression Trees with Model Trees. arXiv preprint arXiv:2006.07493.
}
\author{
Estevão B. Prado: \email{estevao.prado@hotmail.com}.
}

\examples{

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
}
