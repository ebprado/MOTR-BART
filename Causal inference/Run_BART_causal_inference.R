library(grf)
library(bcf)
library(bartBMA)
# devtools::install_github("EoghanONeill/bcfbma")
library(bcfbma)
source('BART_Functions_causal_inference.R')

# data generating process
p = 3 #two control variables and one moderator
n = 250
#
set.seed(1)

x = matrix(rnorm(n*p), nrow=n)

# create targeted selection
q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
# q = rnorm(n)

# generate treatment variable
pi = pnorm(q)
z = rbinom(n,1,pi)

# tau is the true (homogeneous) treatment effect
tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))

# generate the response using q, tau and z
mu = (q + tau*z)
 
# set the noise level relative to the expected mean function of Y
sigma = diff(range(q + tau*pi))/8

# draw the response variable with additive error
y = mu + sigma*rnorm(n)

# If you didn't know pi, you would estimate it here
pihat = pnorm(q)
X_train = as.data.frame(cbind(x,z, pihat))

# My BART implementation for causal inference  -------------------

BART_out = rBART(as.data.frame(cbind(x, pihat)), y, z,
                      num_trees_mu = 5,
                      num_trees_tau = 5,
                      control = list(node_min_size = 5),
                      MCMC = list(iter = 1000,
                                  burn = 1000,
                                  thin = 1))

# bart_causal = apply(BART_out$y_hat,2,mean)
bart_causal = apply(BART_out$tau,2,mean)
plot(tau, bart_causal, main = 'Causal BART'); abline(0,1)
hist(bart_causal, main = 'Causal BART')
mean(bart_causal); sd(bart_causal)

# BART-BMA for causal inference --------- 
pihat = as.matrix(pihat)
bcfBMA_rfunc_example<-bcfBMA(x,y,z,pihat, x.test = x[1:100,], test_pihat = as.matrix(pihat[1:100,]), test_z = z[1:100])
tauhat_bcfbma_linalg<-preds_bcfbma_lin_alg(bcfBMA_rfunc_example,num_iter=2000)
plot(tau[1:100], tauhat_bcfbma_linalg); abline(0,1)
