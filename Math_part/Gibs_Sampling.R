# ----------------------------------------------------------- #
# Description: Gibbs Sampling for Bayesian linear regression  #
#              with conjugate priors both for betas and sigma #
# Author: Estevão Prado                                       #
# Last modification: 22/10/2018                               #
# ----------------------------------------------------------- #
require(MASS)
require(ggplot2)
require(gridExtra)

## MCMC simulation: example 1 (Gibbs sampling)
## -------------------------------------------

## Defining simulation size, burn-in, etc
## --------------------------------------
Niter <- 500
BurnIn <- 2000
TotIter <- Niter+BurnIn
N <- 1000
b0 <- 1.4
b1 <- -0.8
b2 <- 1.5
b3 <- 0.5
betas <- c(b0,b1,b2,b3)
sigma <- 4
AuxBurnIn <- 1

## Simulation scheme 
## -----------------
x1 <- rbinom(n = N, prob = 0.5, size = 1)
x2 <- rpois(n = N, lambda = 5)
x3 <- round(rnorm(n = N, mean = 0, sd = 3),2)
e <- rnorm(n = N, mean = 0, sd = sqrt(sigma))
X <- cbind(rep(1,N),x1, x2, x3)
y <- rnorm(n = N, mean = betas[1] + betas[2]*x1 + betas[3]*x2 + betas[4]*x3 + e)

## Defining prior distributions
## ----------------------------
## Betas ~ Normal(m, V)
## --------------------
m = rbind(0,0,0,0)
V = diag(4)

## Sigma2 ~ Inverse Gamma(a, b)
## ---------------------------
a = 2
b = 1

## Posterior parameters
## --------------------
mu <- solve(t(X)%*%X + solve(V)) %*% (t(X)%*%y + solve(V)%*%m); mu
Lambda <- solve(t(X)%*%X + solve(V)); Lambda
a_ast <- (N/2) + a; a_ast
b_ast <- b + 0.5*(-t(mu)%*%solve(Lambda)%*%mu + t(m)%*%solve(V)%*%m + t(y)%*%y); b_ast

## Data frame that will store MCMC values for betas and sigma2
## -----------------------------------------------------------
SaveResults <- as.data.frame(matrix(data = NA, nrow = Niter, ncol = length(betas)+2))
colnames(SaveResults) <- c('Iter', 'Beta0', 'Beta1', 'Beta2', 'Beta3', 'sigma2')

## Getting started Gibbs Sampling
## ------------------------------
for(i in 1:TotIter){

  MCMCSigma <- 1/rgamma(1, a_ast, b_ast)
  MCMCBetas <- mvrnorm(1, mu, MCMCSigma*Lambda)

  if (i > BurnIn){
    SaveResults[AuxBurnIn,] <- c(AuxBurnIn,MCMCBetas, MCMCSigma)
    AuxBurnIn <- AuxBurnIn + 1
  }
}

## True values
## -----------
betas
sigma

## Plots
## -----
ChainB0 <- ggplot(SaveResults, aes(x=Iter)) +
  geom_line(aes(y=Beta0)) +
  labs(title=expression(paste('Chain of ', beta[0])),
       subtitle='Conjugate priors: Gibbs Sampling',
       x='Iterations',
       y='Values') +
  geom_hline(yintercept=betas[1], linetype = 'dotted', color='coral') +
  theme_bw()

ChainB1 <- ggplot(SaveResults, aes(x=Iter)) +
  geom_line(aes(y=Beta1)) +
  labs(title=expression(paste('Chain of ', beta[1])),
       subtitle='Conjugate priors: Gibbs Sampling',
       x='Iterations',
       y='Values') +
  geom_hline(yintercept=betas[2], linetype = 'dotted', color='coral') +
  theme_bw()

ChainB2 <- ggplot(SaveResults, aes(x=Iter)) +
  geom_line(aes(y=Beta2)) +
  labs(title=expression(paste('Chain of ', beta[2])),
       subtitle='Conjugate priors: Gibbs Sampling',
       x='Iterations',
       y='Values') +
  geom_hline(yintercept=betas[3], linetype = 'dotted', color='coral') +
  theme_bw()

ChainB3 <- ggplot(SaveResults, aes(x=Iter)) +
  geom_line(aes(y=Beta3)) +
  labs(title=expression(paste('Chain of ', beta[3])),
       subtitle='Conjugate priors: Gibbs Sampling',
       x='Iterations',
       y='Values') +
  geom_hline(yintercept=betas[4], linetype = 'dotted', color='coral') +
  theme_bw()

ChainS2 <- ggplot(SaveResults, aes(x=Iter)) +
  geom_line(aes(y=sigma2)) +
  labs(title=expression(paste('Chain of ', sigma^{2})),
       subtitle='Conjugate priors: Gibbs Sampling',
       x='Iterations',
       y='Values') +
  geom_hline(yintercept=sigma, linetype = 'dotted', color='coral') +
  theme_bw()

grid.arrange(ChainB0, ChainB1, ChainB2, ChainB3, ChainS2, ncol=2)