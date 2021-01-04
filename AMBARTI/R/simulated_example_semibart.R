library(devtools)
install_github("ebprado/semibart") # I just added the BART predictions to the output
library(semibart)

install_github("ebprado/MOTR-BART/AMBARTI") # This my implementation of the semibart idea
library(AMBARTI)

################################################################################
# Generate the simulated example
################################################################################

# An Additive Main effects and Multiplicative Interaction (AMMI) effects model

# A Bayesian version of the AMMI model as specified here: https://link.springer.com/content/pdf/10.1007/s13253-014-0168-z.pdf (JosseE et al JABES 2014)
# Andrew Parnell / Danilo Sarti

# In this file, we simulate from the AMMI model specified in the paper above.

rm(list = ls())
library(R2jags)
library(ggplot2)
library(tidyverse)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file

# Likelihood
# Y_{ij} ~ N(mu_{ij}, sigma^2_E)
# with
# mu_{ij} = mu + alpha_i + beta_j + sum_{q=1}^Q lambda_q*gamma_iq*delta_jq
# AA = sum_{q=1}^Q lambda_q*gamma_iq*delta_jq
# Our idea is to estimate alpha_i + beta_j parametrically and the component AA via BART.

# Notation
# Y_ij = response (e.g. yield) for genotype i and environment j, i = 1, ..., I genotypes and j = 1, ..., J environments
# mu is the grand mean
# alpha_i is the genotype effect
# beta_j is the environment effect
# lambda_q is the q-th eigenvalue q = 1,.., Q of the interaction matrix
# Q is the number of components used to model the interaction. Usually Q is fixed at a small number, e.g. 2
# gamma_{iq} is the interaction effect for the q-th eigenvector for genotype i
# delta_{iq} is the interaction effect for the q-th eigenvector for environment j
# E_{ij} is a residual term with E_{ij} ~ N(0, sigma^2_E)
# Usually these models have quite complicated restrictions on the gamma/delta/lambda values but Josse et al show that these are not fully necessary

# Priors
# alpha_i ~ N(0, s_alpha^2)
# beta_j ~ N(0, s_beta^2)

# Simulate data -----------------------------------------------------------

# We will follow the simulation strategy detailed in Section 3.1 of the Josse et al paper

# Specify fixed values
Q = 1 # Number of components
I = 10 # Number of genotypes
J = 10# Number of environments
N = I*J # Total number of obs

# Some further fixed values
mu = 10
sigma_E = 1
alpha = seq(-3,3, length.out = I)
beta = seq(-4,4, length.out = J)
lambda_1 = 12
gamma = seq(2, -2,length.out = I)/sqrt(10)
delta = seq(-0.5, 0.5,length.out = J)

# Now simulate the values
set.seed(123)
G_by_E = expand.grid(1:I, 1:J) ## setting the interaction matrix
mu_ij = mu + alpha[G_by_E[,1]] + beta[G_by_E[,2]]  + lambda_1 * gamma[G_by_E[,1]] * delta[G_by_E[,2]] ## maybe insert lambda2
Y = rnorm(N, mu_ij, sigma_E) ## response variable

##################################
# Semiparametric BART
##################################

# Some pre-processing
x = G_by_E
names(x) = c('g', 'e')
x$g = as.factor(x$g)
x$e = as.factor(x$e)
y = Y

cov_g = x[,'g']
cov_e = x[,'e']

classes_g = sort(unique(cov_g))
classes_e = sort(unique(cov_e))

ng = tapply(cov_g, cov_g, length)
ne = tapply(cov_e, cov_e, length)

x <- model.matrix(~ -1 + g + e, data=x,
                  contrasts.arg=list(g=contrasts(as.factor(x$g), contrasts=F),
                                     e=contrasts(as.factor(x$e), contrasts=F)))
set.seed(001)

# Run Semiparametric BART
semib = semibart(x.train = x, y.train = y, a.train = x)

# Get the main effects estimates
betahat = apply(semib$beta,2,mean)[1:10] # The first 10 are associated to the covariate g (genotype)

# Plot the main effects estimates and add the true values
plot(betahat, cex=2, ylim = c(min(betahat, alpha), max(betahat, alpha)), main='Genotype - semibart') # estimates (black)
points(1:length(alpha), alpha, col=2, cex=2) # true values (red). Looks not too bad.

# Plot the main effects estimates and add the true values
alphahat = apply(semib$beta,2,mean)[11:20] # The remaining 10 are associated to the covariate e (environment)
plot(alphahat, cex=2, ylim = c(min(alphahat, beta), max(alphahat, beta)), main='Environment - semibart') # estimates
points(1:length(beta), beta, col=2, cex=2) # true values. Looks fine.

# Correlation btw y and BART estimate
cor(y, apply(semib$bartfit, 2, mean)) # ~0.31

# Compute the final prediction (y hat)
yhat = x%*%apply(semib$beta,2,mean) + apply(semib$bartfit, 2, mean)
plot(y, yhat, main = 'semibart - y versus y hat'); abline(0,1) # Looks fine
cor(y, yhat); # ~0.89

# Just verifying whether the tree are splitting on the covariates
#aa = t(semib$test) # it correponds to the number of terminal nodes for every tree along the MCMC iteration
#bb = apply(aa, 2, function(x)length(unique(x))) # the number of terminal nodes per tree for each MCMC iteration
#hist(bb, freq = FALSE) # the trees are splitting on the covariates.

##################################
# BART (just to have a benchmark)
##################################
library(BART)
bart = BART::wbart(x, y)
cor(y, bart$yhat.train.mean) # BART and semibart are quite similar. That's fine.

##########################################
# AMBARTI (my implementation of semibart)
##########################################
library(AMBARTI)

# Some pre-processing
x.ambarti = G_by_E
names(x.ambarti) = c('g', 'e')
x.ambarti$g = as.factor(x.ambarti$g)
x.ambarti$e = as.factor(x.ambarti$e)
y = Y
set.seed(101)

# Run AMBARTI (I'm using only 50 trees)
fit.ambarti = ambarti(x.ambarti, y, ntrees = 50, skip_trees = FALSE, nburn = 100, npost = 100, sparse= FALSE)

# Run AMBARTI (I'm using only one tree).
# set.seed(103)
# fit.ambarti = ambarti(x.ambarti, y, ntrees = 1, skip_trees = FALSE, nburn = 100, npost = 100, sparse= FALSE)

# Get the final prediction (y hat)
yhat_ambarti = apply(fit.ambarti$y_hat, 2, mean)
cor(y, yhat_ambarti); # AMBARTI, BART and semibart are quite similar. That's fine.

# Get the prediction specifically from BART
yhat_bart = apply(fit.ambarti$y_hat_bart, 2, mean);
cor(y, yhat_bart); # correlation btw y and BART (AMBARTI package)
cor(y, apply(semib$bartfit, 2, mean)) # correlation btw y and BART estimate (semibart package)

# Plot the main effects estimates and add the true values
plot(1:length(alpha), alpha, col=2, cex=2, main='Genotype - estimates versus true values', ylim= c(-4,4)) # true values
points(apply(fit.ambarti$beta_hat[,1:10], 2, mean), cex=2) # estimates

# Plot the main effects estimates and add the true values
plot(1:length(beta), beta, col=2, cex=2, main='Environment - estimates versus true values') # true values
points(apply(fit.ambarti$beta_hat[,11:20], 2, mean), cex=2, ylim = c(-4,4)) # estimates

# We can see the covariates that are used in the trees have their effect
# in the linear prediction affected (this is easier to be visualise if we
# uncomment rows 153 and 154, and run the results again).

# Tree matrix for the last MCMC iteration. It's possible to see the covariates
# that were used to create the tree structure (column "split_variable").
# From 1:10, the covariates are associated to covariate 'g' (genotypes).
# The remaining ones are associated to covariate 'e' (environment).

fit.ambarti$trees[[100]][[1]] # shows the tree 1 in the last (100) MCMC iteration.

