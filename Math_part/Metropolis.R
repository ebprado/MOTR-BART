# ---------------------------------------------------------------- #
# Description: Metropolis-Hastings for Bayesian linear regression  #
#              with a Laplace distribution as prior for betas.     #
#              In this example, sigma is assumed to be known.      #
# Author: Estevão Prado                                            #
# Last modification: 22/10/2018                                    #
# ---------------------------------------------------------------- #

require(MASS)
require(ggplot2)
require(gridExtra)

## MCMC simulation: example 2 (Metropolis-Hastings)
## ------------------------------------------------

## Defining simulation size, burn-in, etc
## --------------------------------------

Niter <- 10000
BurnIn <- 100
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
x3 <- round(rnorm(n = N, mean = 0, sd = 4),2)
e <- rnorm(n = N, mean = 0, sd = sqrt(sigma))
X <- cbind(rep(1,N),x1, x2, x3)
y <- rnorm(n = N, mean = betas[1] + betas[2]*x1 + betas[3]*x2 + betas[4]*x3 + e)

## Defining prior distributions
## ----------------------------
## Betas ~ Laplace(m_{i}, v_{i})
## ----------------------------
m = rbind(0,0,0,0)
v = rbind(1,1,1,1)

## Data frame that will store MCMC values for betas
## ------------------------------------------------
SaveResults <- as.data.frame(matrix(data = NA, nrow = Niter, ncol = length(betas)+1))
colnames(SaveResults) <- c('Iter', 'Beta0', 'Beta1', 'Beta2', 'Beta3')

## Initial values for Betas and covariance matrix of the proposal distribution
## ---------------------------------------------------------------------------
MCMCBetasI <- c(10,10,10,10)
V = diag(4)*0.0005

## Getting started Metropolis-Hastings
## -----------------------------------
for(i in 1:TotIter){
  
  MCMCBetasC <- mvrnorm(1, MCMCBetasI, V)

  razao <- (-0.5*sigma*(t(y-X%*%MCMCBetasC)%*%(y-X%*%MCMCBetasC) - sum(abs(MCMCBetasC - m)/v))) - 
           (-0.5*sigma*(t(y-X%*%MCMCBetasI)%*%(y-X%*%MCMCBetasI) - sum(abs(MCMCBetasI - m)/v)))
  
  if(runif(1) < min(1, exp(razao)))
    {MCMCBetasI <- MCMCBetasC}
  
  if (i > BurnIn){
    SaveResults[AuxBurnIn,] <- c(AuxBurnIn,MCMCBetasI)
    AuxBurnIn <- AuxBurnIn + 1
  }
}

head(SaveResults)

## True values
## -----------
betas

## Plots
## -----
ChainB0 <- ggplot(SaveResults, aes(x=Iter)) +
  geom_line(aes(y=Beta0)) +
  labs(title=expression(paste('Chain of ', beta[0])),
       subtitle='Non-conjugate priors: Metropolis-Hastings',
       x='Iterations',
       y='Values') +
  geom_hline(yintercept=betas[1], linetype = 'dotted', color='coral') +
  theme_bw()

ChainB1 <- ggplot(SaveResults, aes(x=Iter)) +
  geom_line(aes(y=Beta1)) +
  labs(title=expression(paste('Chain of ', beta[1])),
       subtitle='Non-conjugate priors: Metropolis-Hastings',
       x='Iterations',
       y='Values') +
  geom_hline(yintercept=betas[2], linetype = 'dotted', color='coral') +
  theme_bw()

ChainB2 <- ggplot(SaveResults, aes(x=Iter)) +
  geom_line(aes(y=Beta2)) +
  labs(title=expression(paste('Chain of ', beta[2])),
       subtitle='Non-conjugate priors: Metropolis-Hastings',
       x='Iterations',
       y='Values') +
  geom_hline(yintercept=betas[3], linetype = 'dotted', color='coral') +
  theme_bw()

ChainB3 <- ggplot(SaveResults, aes(x=Iter)) +
  geom_line(aes(y=Beta3)) +
  labs(title=expression(paste('Chain of ', beta[3])),
       subtitle='Non-conjugate priors: Metropolis-Hastings',
       x='Iterations',
       y='Values') +
  geom_hline(yintercept=betas[4], linetype = 'dotted', color='coral') +
  theme_bw()

grid.arrange(ChainB0, ChainB1, ChainB2, ChainB3, ncol=2)
