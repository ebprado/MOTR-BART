#-------------------------------------------------------------------------------------------#
# Description: (BART-Importance Sampling) Generate RMSE for both ATE and CATE based on       #
#              synthetic data described in the simulation section in Hahn et al (2020).     #
#              The results generated here are comparable to those from tables 2 and 3 of    #
#              the aforementioned paper.                                                    #
#-------------------------------------------------------------------------------------------#
options(warn=0)
setwd("~/R/Discussion_paper")
source('Aux_functions.R')

save_file = '~/R/Discussion_paper/results_sim/'
filename = 'IS_Simulation_results'
consolidated_results = NULL
num_rep = 50 # Monte Carlo repetitions
sample = 250
ncov = c(5, 50, 100, 500) # number of covariates
tau_str = c('heterogeneous', 'homogeneous')
mu_str = c('linear', 'nonlinear')

for (s in 1:length(sample)){
  for (k in 1:length(ncov)){
    for (m in 1:length(mu_str)){
      for (t in 1:length(tau_str)){
        for (i in 1:num_rep){

          print(paste('tau_str = ', tau_str[t]))
          print(paste('mu_str = ', mu_str[m]))
          print(paste('p = ', ncov[k]))
          print(paste('n = ', sample[s]))
          print(paste('rep = ', i))

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

          # BCF-IS -----------------------------------------------------------------------------------------

          j = 0
          fit.safeBCF = NULL
          class(fit.safeBCF) = NULL

          while(j <= 3 & any(class(fit.safeBCF) %in% c('error', 'NULL'))) {
            fit.safeBCF = tryCatch({
              beta_par <- 1


              lambda_mu <- 0.45
              lambda_tau <- 0.05


              Num_models <- 3000
              num_trees1_mu <- 50
              num_trees1_tau <- 25

              seed1 <- 42
              ncores <- 4


              safeBCF_with_intervals(seed=seed1,
                                     y=y,
                                     original_datamat=x,
                                     ztrain=z,
                                     pihat_train=as.matrix(pihat),
                                     #test_datamat=NULL,
                                     #test_z=NULL,
                                     #test_pihat=NULL,
                                     lambda_mu=lambda_mu,
                                     lambda_tau=lambda_tau,
                                     num_models=Num_models,
                                     num_trees_mu=num_trees1_mu,
                                     num_trees_tau=num_trees1_tau,
                                     beta_par=beta_par,
                                     ncores=ncores,
                                     outsamppreds=1,
                                     nu=3,
                                     a_mu=5,
                                     a_tau=5,
                                     sigquant=0.9,
                                     valid_trees=1,
                                     tree_prior=1,
                                     imp_sampler=1,
                                     alpha_BCF_mu=0.95,
                                     beta_BCF_mu=2,
                                     alpha_BCF_tau=0.25,
                                     beta_BCF_tau=3,
                                     include_pi= "control",
                                     fast_approx = 0,
                                     PIT_propensity=0,
                                     l_quant=0.025,
                                     u_quant=0.975,
                                     root_alg_precision=0.00001)


            },error = function(e) e)

            j = j+1
          }

          if (class(fit.safeBCF)[1] != 'simpleError'){


            beta_par <- 1


            lambda_mu <- 0.45
            lambda_tau <- 0.05


            Num_models <- 3000
            num_trees1_mu <- 50
            num_trees1_tau <- 25

            seed1 <- 42
            ncores <- 4


            testsafeBCF <- safeBCF_with_intervals(seed=seed1,
                                                  y=y,
                                                  original_datamat=x,
                                                  ztrain=z,
                                                  pihat_train=as.matrix(pihat),
                                                  test_datamat=t.x,
                                                  # test_z=t.z,
                                                  test_pihat=as.matrix(t.pihat),
                                                  lambda_mu=lambda_mu,
                                                  lambda_tau=lambda_tau,
                                                  num_models=Num_models,
                                                  num_trees_mu=num_trees1_mu,
                                                  num_trees_tau=num_trees1_tau,
                                                  beta_par=beta_par,
                                                  ncores=ncores,
                                                  outsamppreds=1,
                                                  nu=3,
                                                  a_mu=5,
                                                  a_tau=5,
                                                  sigquant=0.9,
                                                  valid_trees=1,
                                                  tree_prior=1,
                                                  imp_sampler=1,
                                                  alpha_BCF_mu=0.95,
                                                  beta_BCF_mu=2,
                                                  alpha_BCF_tau=0.25,
                                                  beta_BCF_tau=3,
                                                  include_pi= "control",
                                                  fast_approx = 0,
                                                  PIT_propensity=0,
                                                  l_quant=0.025,
                                                  u_quant=0.975,
                                                  root_alg_precision=0.00001)

            tau.hat.safeBCF.test = fit.safeBCF$ITEests
            aux_safeBCF = c('BCF-IS',
                            rmse_cate(tau, fit.safeBCF$ITEests),
                            rmse_ate(tau, fit.safeBCF$CATEests),
                            rmse_cate(t.tau, testsafeBCF$ITEests),
                            rmse_ate(t.tau, testsafeBCF$CATEests))

          } else {
            aux_safeBCF = c('BCF-IS',
                            NA,
                            NA,
                            NA,
                            NA)
          }

          # BART-IS -----------------------------------------------------------------------------------------

          j = 0
          fit.safeBART = NULL
          class(fit.safeBART) = NULL

          while(j <= 3 & any(class(fit.safeBART) %in% c('error', 'NULL'))) {
            fit.safeBART = tryCatch({
              beta_par <- 1


              lambda_mu <- 0.45
              lambda_tau <- 0.05


              Num_models <- 3000
              num_trees2 <- 30

              seed1 <- 42
              ncores <- 4


              safeBart_ITEs(seed=seed1,
                            y=y,
                            original_datamat=x,
                            ztrain=z,
                            pihat_train=as.matrix(pihat),
                            #test_datamat=NULL,
                            #test_z=NULL,
                            #test_pihat=NULL,
                            num_models=Num_models,
                            num_trees=num_trees2,
                            beta_par=beta_par,
                            ncores=ncores,
                            outsamppreds=1,
                            nu=3,
                            a=10,
                            sigquant=0.9,
                            valid_trees=1,
                            tree_prior=1,
                            imp_sampler=1,
                            alpha_BART=0.95,
                            beta_BART=2,
                            fast_approx = 0,
                            PIT_propensity=0,
                            l_quant=0.025,
                            u_quant=0.975,
                            root_alg_precision=0.00001)

            },error = function(e) e)

            j = j+1
          }

          if (class(fit.safeBART)[1] != 'simpleError'){


            beta_par <- 1


            lambda_mu <- 0.45
            lambda_tau <- 0.05


            Num_models <- 3000
            num_trees2 <- 30

            seed1 <- 42
            ncores <- 4


            testsafeBART <- safeBart_ITEs(seed=seed1,
                                          y=y,
                                          original_datamat=x,
                                          ztrain=z,
                                          pihat_train=as.matrix(pihat),
                                          test_datamat=t.x,
                                          # test_z=t.z,
                                          test_pihat=as.matrix(t.pihat),
                                          num_models=Num_models,
                                          num_trees=num_trees2,
                                          beta_par=beta_par,
                                          ncores=ncores,
                                          outsamppreds=1,
                                          nu=3,
                                          a=10,
                                          sigquant=0.9,
                                          valid_trees=1,
                                          tree_prior=1,
                                          imp_sampler=1,
                                          alpha_BART=0.95,
                                          beta_BART=2,
                                          fast_approx = 0,
                                          PIT_propensity=0,
                                          l_quant=0.025,
                                          u_quant=0.975,
                                          root_alg_precision=0.00001)

            tau.hat.safeBART.test = fit.safeBART$ITEests
            aux_safeBART = c('BART-IS',
                             rmse_cate(tau, fit.safeBART$ITEests),
                             rmse_ate(tau, fit.safeBART$CATEest),
                             rmse_cate(t.tau, testsafeBART$ITEests),
                             rmse_ate(t.tau, testsafeBART$CATEest))

          } else {
            aux_safeBART = c('BART-IS',
                             NA,
                             NA,
                             NA,
                             NA)
          }


          # Save results

          SaveResults = as.data.frame(
            rbind(
              aux_safeBCF,
              aux_safeBART
            )
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