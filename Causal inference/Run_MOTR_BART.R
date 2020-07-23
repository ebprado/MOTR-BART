#-------------------------------------------------------------------------------------------#
# Description: (MOTR-BART) Generate RMSE for both ATE and CATE based on synthetic data       #
#              described in the simulation section in Hahn et al (2020). The results        #
#              generated here are comparable to those from tables 2 and 3 of the            #
#              aforementioned paper.                                                        #
#-------------------------------------------------------------------------------------------#

options(warn=0)
setwd("~/R/Discussion_paper")
source('Aux_functions.R')
source('MOTR_BART.R')

save_file = '~/R/Discussion_paper/results_sim/'
filename = '250MOTR_BART_Simulation_results'
consolidated_results = NULL
num_rep = 50 # Monte Carlo repetitions
sample = 250 # sample
ncov = c(5, 50, 100, 500) # number of covariates
tau_str = c('heterogeneous', 'homogeneous')
mu_str = c('linear', 'nonlinear')

for (s in 1:length(sample)){
  for (k in 1:length(ncov)){
    for (m in 1:length(mu_str)){
      for (t in 1:length(tau_str)){
        for (i in 1:num_rep){

          # Train data -----------
          data = generate_data(n = sample[s], p = ncov[k], tau = tau_str[t], mu = mu_str[m], seed = i + k)

          p = data$p
          n = data$n
          y = data$y
          x = data$x
          z = data$z
          pihat = data$pihat
          tau = data$tau

          x.mod = data.frame(x,z,pihat)
          x.train = makeModelMatrixFromDataFrame(data.frame(rbind(x,x),z = c(rep(1,n),rep(0,n)), pihat = c(pihat, pihat)))

          # Test data -----------------

          data.test = generate_data(n = sample[s], p = ncov[k], tau = tau_str[t], mu = mu_str[m], seed = 1000 + i + k)
          t.y = data.test$y
          t.x = data.test$x
          t.z = data.test$z
          t.pihat = data.test$pihat
          t.tau = data.test$tau

          t.x.mod = data.frame(t.x,t.z,t.pihat)
          x.test = makeModelMatrixFromDataFrame(data.frame(rbind(t.x,t.x), z = c(rep(1,n),rep(0,n)), pihat = c(t.pihat, t.pihat)))

          # MOTR-BART -----------------------------------------------------------------------------------------
          motr_bart = rBART(x.mod, y, num_trees = 10,
                            control = list(node_min_size = 5),
                            MCMC = list(iter = 1000,
                                        burn = 4000,
                                        thin = 1))

          aux.tau.motr = predict_rBART(x.train, motr_bart, type = 'mean')
          tau.motr.bart = aux.tau.motr[1:n] -  aux.tau.motr[1:n+n]
          test.aux.tau.motr = predict_rBART(x.test, motr_bart, type = 'mean')
          test.tau.motr = test.aux.tau.motr[1:n] -  test.aux.tau.motr[1:n+n]
          aux_motrbart = c('ps-MOTR-BART',
                           rmse_cate(tau, tau.motr.bart),
                           rmse_ate(tau, tau.motr.bart),
                           rmse_cate(t.tau, test.tau.motr),
                           rmse_ate(t.tau, test.tau.motr))

          # Save results

          SaveResults = as.data.frame(
            rbind(aux_motrbart)
          )

          SaveResults$rep = i
          SaveResults$n = n
          SaveResults$p = p
          SaveResults$tau_str = tau_str[t]
          SaveResults$mu_str = mu_str[m]

          consolidated_results = rbind(consolidated_results, SaveResults)
          # rownames(consolidated_results) = NULL
          # names(consolidated_results) = c('Algorithm', 'RMSE_CATE_train', 'RMSE_ATE_train', 'RMSE_CATE_test', 'RMSE_ATE_test', 'rep', 'n','tau_str', 'mu_str')
          save(consolidated_results, file=paste(save_file, filename,'.RData', sep=''))
        }
      }
    }
  }
}