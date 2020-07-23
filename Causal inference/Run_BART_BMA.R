#-------------------------------------------------------------------------------------------#
# Description: (BART-BMA) Generate RMSE for both ATE and CATE based on                      #
#              synthetic data described in the simulation section in Hahn et al (2020).     #
#              The results generated here are comparable to those from tables 2 and 3 of    #
#              the aforementioned paper.                                                    #
#-------------------------------------------------------------------------------------------#
options(warn=0)
setwd("~/R/Discussion_paper")
source('Aux_functions.R')
# source('MOTR_BART.R')

save_file = '~/R/Discussion_paper/results_sim/'
filename = 'BMA_Simulation_results'
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

# BART-BMA-----------------------------------------------------------------------------------------

j = 0
fit.bartBMA = NULL
class(fit.bartBMA) = NULL

while(j <= 3 & any(class(fit.bartBMA) %in% c('error', 'NULL'))) {
  fit.bartBMA = tryCatch({
    #bartBMA(x.mod, y, x.test = x.train)

    temp_x_pihat <- cbind(x,pihat)
    bartBMA_with_ITEs_exact_par(0.025,0.975,newdata=NULL,update_resids=1,
                                num_cores=1,root_alg_precision=0.00001,
                                temp_x_pihat,z ,y,
                                a=3,nu=3,sigquant=0.9,c=1000,
                                pen=12,num_cp=20,#x.test=matrix(0.0,0,0),
                                num_rounds=5,alpha=0.95,beta=2,
                                split_rule_node=0,
                                gridpoint=1,maxOWsize=100,
                                num_splits=5,gridsize=10,zero_split=1,
                                only_max_num_trees=1,
                                min_num_obs_after_split = 5,
                                exact_residuals = 1)




  },error = function(e) e)

  j = j+1
}

if (class(fit.bartBMA)[1] != 'simpleError'){

  tau.hat.bartBMA = fit.bartBMA[[2]]
  #test.bartBMA.prediction = predict_bartBMA(fit.bartBMA, x.test)

  test_bbma_temp <- bartBMA_with_ITEs_exact_par(0.025,0.975,
                                                newdata=cbind(t.x,t.pihat),update_resids=1,
                                                num_cores=1,root_alg_precision=0.00001,
                                                cbind(x,pihat),z ,y,
                                                a=3,nu=3,sigquant=0.9,c=1000,
                                                pen=12,num_cp=20,#x.test1=matrix(0.0,0,0),
                                                num_rounds=5,alpha=0.95,beta=2,
                                                split_rule_node=0,
                                                gridpoint=1,maxOWsize=100,
                                                num_splits=5,gridsize=10,zero_split=1,
                                                only_max_num_trees=1,
                                                min_num_obs_after_split = 5,
                                                exact_residuals = 1)

  tau.hat.bartBMA.test = test_bbma_temp[[2]]
  aux_bartBMA = c('ps-BART-BMA',
                  rmse_cate(tau, tau.hat.bartBMA),
                  sqrt((mean(tau) - fit.bartBMA[[3]])^2),
                  rmse_cate(t.tau, tau.hat.bartBMA.test),
                  sqrt((mean(t.tau) - test_bbma_temp[[3]])^2))

} else {
  aux_bartBMA = c('ps-BART-BMA',
                  NA,
                  NA,
                  NA,
                  NA)
}

# Save results

SaveResults = as.data.frame(
  rbind(
        aux_bartBMA
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