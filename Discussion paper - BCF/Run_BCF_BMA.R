#-------------------------------------------------------------------------------------------#
# Description: (BCF-BMA) Generate RMSE for both ATE and CATE based on                       #
#              synthetic data described in the simulation section in Hahn et al (2020).     #
#              The results generated here are comparable to those from tables 2 and 3 of    #
#              the aforementioned paper.                                                    #
#-------------------------------------------------------------------------------------------#
options(warn=0)
setwd("~/R/Discussion_paper")
source('Aux_functions.R')
# source('MOTR_BART.R')

save_file = '~/R/Discussion_paper/results_sim/'
filename = 'BCF_BMA_Simulation_results'
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

          # BCF-BMA-----------------------------------------------------------------------------------------

          j = 0
          fit.bcfBMA = NULL
          class(fit.bcfBMA) = NULL

          while(j <= 3 & any(class(fit.bcfBMA) %in% c('error', 'NULL'))) {
            fit.bcfBMA = tryCatch({
              bcfBMA(x,y,z,pihat,
                     a_mu=1,a_tau=0.5,nu=3,sigquant=0.9,c=100,
                     pen_mu=12,pen_tau=12,num_cp_mu=20,num_cp_tau=20,
                     # x.test=t.x,test_z =t.z,
                     # test_pihat = as.matrix(t.pihat),
                     ntree_control=5,ntree_moderate=5,
                     alpha_mu=0.95,alpha_tau=0.25,beta_mu=2,beta_tau=3,split_rule_node=1,
                     gridpoint=1,maxOWsize=100,num_splits_mu =5, num_splits_tau =5,
                     gridsize_mu=20, gridsize_tau=20, include_pi= "control",
                     zero_split = 1, only_max_num_trees = 1,mu_or_tau_each_round = 1,separate_tree_numbers = 1,
                     min_num_obs_after_mu_split = 5, min_num_obs_after_tau_split = 5,
                     transform_resids = 0)

            },error = function(e) e)

            j = j+1
          }



          j = 0
          fit.bcfBMA22 = NULL
          class(fit.bcfBMA22) = NULL

          while(j <= 3 & any(class(fit.bcfBMA22) %in% c('error', 'NULL'))) {
            fit.bcfBMA22 = tryCatch({
              bcfBMA(x,y,z,pihat,
                     a_mu=1,a_tau=0.5,nu=3,sigquant=0.9,c=100,
                     pen_mu=12,pen_tau=12,num_cp_mu=20,num_cp_tau=20,
                     x.test=t.x,test_z =t.z,
                     test_pihat = as.matrix(t.pihat),
                     ntree_control=5,ntree_moderate=5,
                     alpha_mu=0.95,alpha_tau=0.25,beta_mu=2,beta_tau=3,split_rule_node=1,
                     gridpoint=1,maxOWsize=100,num_splits_mu =5, num_splits_tau =5,
                     gridsize_mu=20, gridsize_tau=20, include_pi= "control",
                     zero_split = 1, only_max_num_trees = 1,mu_or_tau_each_round = 1,separate_tree_numbers = 1,
                     min_num_obs_after_mu_split = 5, min_num_obs_after_tau_split = 5,
                     transform_resids = 0)

            },error = function(e) e)

            j = j+1
          }



          if (class(fit.bcfBMA)[1] != 'simpleError' & class(fit.bcfBMA22)[1] != 'simpleError'){
            #tau.hat.bcfBMA <- preds_bcfbma_lin_alg(fit.bcfBMA, num_iter=2000)
            #tau.hat.bcfBMA <- fit.bcfBMA$ITE_estimates
            bcfBMAtrain <- pred_ints_exact_bcf_TE(fit.bcfBMA,l_quant=0.025,u_quant=0.975,num_cores=1,root_alg_precision = 0.00001)
            bcfBMAtest <- pred_ints_exact_bcf_TE(fit.bcfBMA22,l_quant=0.025,u_quant=0.975,newdata=t.x,num_cores=1,root_alg_precision = 0.00001)
            #test.bcfBMA.prediction <- fit.bcfBMA$test.preds_tau

            aux_bcfBMA = c('BCF-BMA',
                           rmse_cate(tau, bcfBMAtrain$ITE_estimates),
                           rmse_ate(tau, bcfBMAtrain$CATE_estimate),
                           rmse_cate(t.tau, bcfBMAtest$ITE_estimates),
                           rmse_ate(t.tau, bcfBMAtest$CATE_estimate))

          } else {
            aux_bcfBMA = c('BCF-BMA',
                           NA,
                           NA,
                           NA,
                           NA)
          }

          # Save results

          SaveResults = as.data.frame(
            rbind(
              aux_bcfBMA
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