library(dplyr)
library(ggplot2)
library(xtable)
library(tidyverse)
library(stringr)
options(warn = 0)

setwd("~/R/Discussion_paper/results_sim")

# Load the results ----------------------------------------------------------------------------------------

load('BART_CF_Simulation_results.RData'); res1 = consolidated_results
load('BCF_Simulation_results.RData'); res2 = consolidated_results
load('1BMA_Simulation_results.RData'); res3 = consolidated_results
load('2BMA_Simulation_results.RData'); res4 = consolidated_results
load('1250MOTR_BART_Simulation_results.RData'); res5 = consolidated_results
load('2250MOTR_BART_Simulation_results.RData'); res6 = consolidated_results
load('3250MOTR_BART_Simulation_results.RData'); res7 = consolidated_results
load('4250MOTR_BART_Simulation_results.RData'); res8 = consolidated_results
# load('1BCF_BMA_Simulation_results.RData'); res9 = consolidated_results
# load('2BCF_BMA_Simulation_results.RData'); res10 = consolidated_results
db = as.data.frame(rbind(res1, res2, res3, res4, res5, res6, res7, res8))

# Organise some details ----------------------------------------------------------------------------------------

rownames(db) = NULL
names(db) = c('Algorithm', 'RMSE_CATE_train', 'RMSE_ATE_train', 'RMSE_CATE_test', 'RMSE_ATE_test', 'rep', 'n', 'p', 'tau_str', 'mu_str')
db$RMSE_CATE_train = as.numeric(as.character(db$RMSE_CATE_train))
db$RMSE_ATE_train = as.numeric(as.character(db$RMSE_ATE_train))
db$RMSE_CATE_test = as.numeric(as.character(db$RMSE_CATE_test))
db$RMSE_ATE_test = as.numeric(as.character(db$RMSE_ATE_test))
db$rep = as.numeric(as.character(db$rep))
db$p = factor(db$p, levels = c('5', '50', '100', '500'))
db$tau_str = factor(db$tau_str, levels = c('heterogeneous', 'homogeneous'), labels = c(expression(paste(tau, '(x) =  heterogeneous')), expression(paste(tau, '(x) = homogeneous'))))
db$mu_str = factor(db$mu_str, levels = c('linear', 'nonlinear'), labels = c(expression(paste(mu, '(x) = linear')), expression(paste(mu, '(x) = nonlinear'))))
db$Algorithm = factor(db$Algorithm, levels=c('Causal Forests','BCF', 'ps-BART-BMA', 'ps-DART', 'ps-BART', 'ps-MOTR-BART'), labels=c('Causal RF', 'BCF', 'ps-BART-BMA', 'ps-DART', 'ps-BART', 'ps-MOTR BART'))
# db$Algorithm = factor(db$Algorithm, levels=c('Causal Forests','BCF', 'BCF-BMA', 'ps-BART-BMA', 'ps-DART', 'ps-BART', 'ps-MOTR-BART'), labels=c('Causal RF', 'BCF', 'BCF-BMA', 'ps-BART-BMA', 'ps-DART', 'ps-BART', 'ps-MOTR BART'))

# Plot results ------------------------------------------------------------------------------------------
plott = function(sample_size, treat_effect, estimand, order){

  pdf(paste(estimand, '.pdf', sep = ''), width = 10, height = 8)

  my_plot <- db %>%
    dplyr::filter(n == sample_size) %>%
    ggplot(aes_string(x = 'p', y=paste('RMSE_', estimand, sep=''), colour='Algorithm')) +
    geom_boxplot() +
    scale_color_discrete("Algorithm", drop=FALSE) +
    labs(title='',
         y = 'RMSE') +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          legend.position = 'bottom',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          strip.background = element_rect(fill="white", size = 0.8)
    )  +
    facet_grid(mu_str ~ tau_str, scales = 'free', labeller = label_parsed) +
    labs(colour='')

  print(my_plot)
  dev.off()

}

plott(250, 'homogeneous', 'CATE_train', '1') # 10 x 8
plott(250, 'homogeneous', 'ATE_train', '5') # 10 x 8
# plott(250, 'homogeneous', 'CATE_test', '1') # 10 x 8
# plott(250, 'homogeneous', 'ATE_test', '5') # 10 x 8
