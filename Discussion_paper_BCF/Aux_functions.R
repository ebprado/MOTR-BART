#------------------------------------------------------------------------------#
# Author: Estevao Prado                                                        #
# Description: Functions to compute RMSE for ATE and CATE as well as to        #
#              to generate the data following Hahn et al (2020)                #
# Last modified date: 23/07/2020                                               #
#------------------------------------------------------------------------------#
library(dbarts)
library(bcf)
library(grf)
library(bayeslm)
library(nnet)
library(bartBMA)
library(bcfbma)
library(safeBart)
library(BART)

g = function(x){
  ifelse(x == 1, 2, ifelse(x == 2, -1, -4))
}

rmse_ate = function(tau, tau_hat){
  return(sqrt((mean(tau) - mean(tau_hat))^2))
}

rmse_cate = function(tau, tau_hat){
  return(sqrt(mean((tau - tau_hat)^2)))
}

# Simulating data sets following Hahn et al (2020)
generate_data = function(n, p, tau, mu, seed){

  set.seed(seed)

  x = matrix(rnorm(n*p), nrow=n, ncol=p)
  x[,1] = rnorm(n)
  x[,2] = rnorm(n)
  x[,3] = rnorm(n)
  x[,4] = rbinom(n, 1, 0.5)
  x[,5] = sample(c(1,2,3), n, replace = TRUE)
  u = runif(n,0,1)

  # Tau --------
  if (tau == 'homogeneous'){tau_x = 3}
  if (tau == 'heterogeneous'){tau_x = 1 + 2*x[,2]*x[,4]} # x[,4] after talking to Eoghan. It makes sense!

  # Mu -----------
  if (mu == 'linear'){mu_x = 1 + g(x[,5]) + x[,1]*x[,3]}
  if (mu == 'nonlinear'){ mu_x = -6 + g(x[,5]) + 6*abs(x[,3] - 1)}

  # Pi and Pihat --------
  pi.x = 0.8*pnorm((3*mu_x/sd(mu_x)) - 0.5*x[,1]) + 0.05 + u/10
  z = rbinom(n, 1, pi.x)
  pihat = apply(bart(x,z, verbose=FALSE)$yhat.train, 2, mean)

  # Response variable
  y = mu_x + tau_x * z

  return(list(y = y,
              x = x,
              z = z,
              tau = tau_x,
              mu = mu_x,
              pihat = pihat,
              n = n,
              p = p))
}