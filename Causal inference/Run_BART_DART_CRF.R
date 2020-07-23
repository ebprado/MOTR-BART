#-------------------------------------------------------------------------------------------#
# Description: (BART, DART and Causal RF) Generate RMSE for both ATE and CATE based on      #
#              synthetic data described in the simulation section in Hahn et al (2020).     #
#              The results generated here are comparable to those from tables 2 and 3 of    #
#              the aforementioned paper.                                                    #
#-------------------------------------------------------------------------------------------#
options(warn=0)
setwd("~/R/Discussion_paper")
source('Aux_functions.R')
source('MOTR_BART.R')

save_file = '~/R/Discussion_paper/results_sim/'
filename = 'BART_CF_Simulation_results'
consolidated_results = NULL
num_rep = 50 # Monte Carlo repetitions
sample = c(250) # sample
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

          # BART-----------------------------------------------------------------------------------------

          j = 0
          fit.bart = NULL
          class(fit.bart) = NULL

          while(j <= 3 & any(class(fit.bart) %in% c('error', 'NULL'))) {

            fit.bart = tryCatch({
              wbart(x.mod, y, x.train, ndpost=10000L, nskip = 1000, sigest = sd(y))
            },error = function(e) e)

            j = j+1
          }

          if (class(fit.bart)[1] != 'simpleError'){

            tau.hat.bart = fit.bart$yhat.test.mean[1:n] -  fit.bart$yhat.test.mean[1:n+n]
            test.bart.prediction = apply(predict(fit.bart, x.test), 2, mean)
            tau.hat.bart.test = test.bart.prediction[1:n] - test.bart.prediction[1:n+n]
            aux_bart = c('ps-BART',
                         rmse_cate(tau, tau.hat.bart),
                         rmse_ate(tau, tau.hat.bart),
                         rmse_cate(t.tau, tau.hat.bart.test),
                         rmse_ate(t.tau, tau.hat.bart))

          } else {
            aux_bart = c('ps-BART',
                         NA,
                         NA,
                         NA,
                         NA)
          }

          # Dirichlet - BART (DART)-----------------------------------------------------------------------------------------

          j = 0
          fit.dart = NULL
          class(fit.dart) = NULL

          while(j <= 3 & any(class(fit.dart) %in% c('error', 'NULL'))) {

            fit.dart = tryCatch({
              wbart(x.mod, y, x.train, sparse = TRUE, ndpost=10000L, nskip = 1000, sigest = sd(y))
            },error = function(e) e)

            j = j+1
          }

          if (class(fit.bart)[1] != 'simpleError'){

            tau.hat.dart = fit.dart$yhat.test.mean[1:n] -  fit.dart$yhat.test.mean[1:n+n]
            test.dart.prediction = apply(predict(fit.dart, x.test), 2, mean)
            tau.hat.dart.test = test.dart.prediction[1:n] - test.dart.prediction[1:n+n]
            aux_dart = c('ps-DART',
                         rmse_cate(tau, tau.hat.dart),
                         rmse_ate(tau, tau.hat.dart),
                         rmse_cate(t.tau, tau.hat.dart.test),
                         rmse_ate(t.tau, tau.hat.dart))

          } else {
            aux_dart = c('ps-DART',
                         NA,
                         NA,
                         NA,
                         NA)
          }

          # Causal Forest -----------------------------------------------------------------------------------------
          fit.causalForest = causal_forest(x, y, z, num.trees = 4000)
          tau.hat.causalForest = fit.causalForest$predictions
          test = average_treatment_effect(fit.causalForest, target.sample = 'all')
          test.causalForest.prediction = as.matrix(predict(fit.causalForest, newdata = t.x))
          aux_causalForest = c('Causal Forests',
                               rmse_cate(tau, tau.hat.causalForest),
                               rmse_ate(tau, tau.hat.causalForest),
                               rmse_cate(t.tau, test.causalForest.prediction),
                               rmse_ate(t.tau, test.causalForest.prediction))


          # Save results

          SaveResults = as.data.frame(
            rbind(aux_bart,
                  aux_dart,
                  aux_causalForest)
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