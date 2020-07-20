library(dplyr)
library(ggplot2)
library(xtable)
library(tidyverse)
library(stringr)
options(warn = 0)
# setwd("~/R/Discussion_paper/previous_results")
setwd("~/R/Discussion_paper/results_sim")
load('BART_CF_Simulation_results.RData'); res1 = consolidated_results
load('BCF_Simulation_results.RData'); res2 = consolidated_results
load('1BMA_Simulation_results.RData'); res3 = consolidated_results
load('2BMA_Simulation_results.RData'); res4 = consolidated_results
load('1250MOTR_BART_Simulation_results.RData'); res5 = consolidated_results
load('2250MOTR_BART_Simulation_results.RData'); res6 = consolidated_results
load('3250MOTR_BART_Simulation_results.RData'); res7 = consolidated_results
load('4250MOTR_BART_Simulation_results.RData'); res8 = consolidated_results
load('1BCF_BMA_Simulation_results.RData'); res9 = consolidated_results
load('2BCF_BMA_Simulation_results.RData'); res10 = consolidated_results
# load('IS_Simulation_results.RData'); res11 = consolidated_results

db = as.data.frame(rbind(res1, res2, res3, res4, res5, res6, res7, res8, res9, res10))
# db = as.data.frame(rbind(res3, res4))
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
db$Algorithm = factor(db$Algorithm, levels=c('Causal Forests','BCF', 'BCF-BMA', 'ps-BART-BMA', 'ps-DART', 'ps-BART', 'ps-MOTR-BART'), labels=c('Causal RF', 'BCF', 'BCF-BMA', 'ps-BART-BMA', 'ps-DART', 'ps-BART', 'ps-MOTR BART'))
# db$Algorithm = factor(db$Algorithm, levels=c('Causal Forests','BCF', 'BCF-BMA', 'BCF-IS', 'BART-IS', 'ps-BART-BMA', 'ps-DART', 'ps-BART', 'ps-MOTR-BART'), labels=c('Causal RF', 'BCF', 'BCF-BMA', 'BCF-IS', 'BART-IS', 'ps-BART-BMA', 'ps-DART', 'ps-BART', 'ps-MOTR BART'))

plott = function(sample_size, treat_effect, estimand, order){
  
pdf(paste(paste0(order, estimand, '_', treat_effect, '_', sample_size,'.pdf', sep = ''), sep = ''), width = 6, height = 10)
  
my_plot <- db %>%
  dplyr::filter(n == sample_size, str_detect(tau_str, treat_effect)) %>% 
  ggplot(aes_string(x = 'p', y=paste('RMSE_', estimand, sep=''), colour='Algorithm')) + 
  geom_boxplot() +
  scale_color_discrete("Algorithm", drop=FALSE, ) +
  labs(#title = '(a)',
    title = paste(estimand, ' tau(x) = ',treat_effect),
    # subtitle = 'Descriptive statistics',
    # x = 'Dataset',
    y = 'RMSE') +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        legend.position = 'bottom',
        # axis.text.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank()
        # panel.background = element_blank()
  )  +
  # facet_grid(. ~ num_trees_class, scales='free')
  # facet_grid(Dataset ~ ., scales='free')
  facet_wrap(. ~ mu_str , scales='free', ncol=1, labeller = label_parsed) +
  # facet_grid(mu_str ~ tau_str, scales = 'free') +
  # facet_wrap( ~ tau_str + mu_str, scales='free') +
  labs(colour='')

print(my_plot)
dev.off()

}

plott(250, 'homogeneous', 'CATE_train', '1')
plott(250, 'heterogeneous', 'CATE_train', '3')
plott(250, 'homogeneous', 'ATE_train', '5')
plott(250, 'heterogeneous', 'ATE_train', '7')

plott(250, 'homogeneous', 'CATE_test', '2')
plott(250, 'heterogeneous', 'CATE_test', '4')
plott(250, 'homogeneous', 'ATE_test', '6')
plott(250, 'heterogeneous', 'ATE_test', '8')
