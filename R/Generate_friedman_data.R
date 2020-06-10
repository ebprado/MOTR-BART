# Simple version ----------------------------------------------------------

sim_friedman_simple = function(n, p = 0, scale_err = 1) {
  # Simulate some data using a multivariate version of Friedman
  # y = 10sin(πx1x2)+20(x3−0.5)2+10x4+5x5+ε
  
  X = matrix(runif(n*(5+p)), nrow = n, ncol = 5 + p)
  
  pars = c(10, 20, 10, 5)
  
  mean = pars[1] * sin(pi*X[,1]*X[,2]) + pars[2] * (X[,3]-0.5)^2 + 
    pars[3] * X[,4] + pars[4] * X[,5]
  
  y = rnorm(n, mean, scale_err)
  
  file = paste('~/R/simulated_data/',
               paste('friedman','n' , n, 'p', p+5, sep='_'),'.txt', sep='')
  
  write.table(cbind(X, y), file=file)
}

set.seed(001)
sim_friedman_simple(200, 0)
sim_friedman_simple(200, 5)
sim_friedman_simple(200, 45)
sim_friedman_simple(500, 0)
sim_friedman_simple(500, 5)
sim_friedman_simple(500, 45)
sim_friedman_simple(1000, 0)
sim_friedman_simple(1000, 5)
sim_friedman_simple(1000, 45)
