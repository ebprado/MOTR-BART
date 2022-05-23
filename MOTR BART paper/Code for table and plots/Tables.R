library(xtable)
library(tidyr)
library(reshape2)
####################################################################################
## Friedman dataset
####################################################################################

setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_Friedman_DBs/Simulation_results/")
# setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_Friedman_DBs_10_trees_1000_burn_5000_post/Simulation_results/")

load('Results_Friedman_DBs.RData')

tab = save_results_friedman
tab$Algorithm = factor(tab$Algorithm, levels=c('CVGLMnet', 'RF', 'GB', 'BART_def', 'BART', 'LM_BART'), labels=c('Lasso', 'RF','GB', 'BART (200 trees)', 'BART (10 trees)', 'MOTR BART'))
tab$obs = 1:(dim(tab))[1]
tab$sample_size = factor(sprintf("%02d", as.numeric(str_split_fixed(tab$Dataset, '_', n = 5)[,3])), levels=c('200', '500', '1000'), labels=c('n = 200', 'n = 500', 'n = 1000'))
tab$num_cov = factor(sprintf("%02d", as.numeric(str_split_fixed(tab$Dataset, '_', n = 5)[,5])), level=c('05', '10', '50'), labels=c('5', '10', '50'))
head(tab)
tab$num_trees_class <- factor(tab$num_trees, levels = c('10', '200', '<NA>'), labels = c('T = 10', 'T = 200', 'No trees'))

################################################################################
## Median (1st and 3rd quartiles) of the train and test data
################################################################################

table1 <- tab %>%
  group_by(Algorithm, num_cov, sample_size) %>%
  summarise(RMSE_train_median = round(median(RMSE_train), 2),
            RMSE_train_1st = sprintf("%.2f", round(quantile(RMSE_train, 0.25), 2)),
            RMSE_train_3rd = sprintf("%.2f", round(quantile(RMSE_train, 0.75), 2)),
            RMSE_test_median = round(median(RMSE_test), 2),
            RMSE_test_1st = sprintf("%.2f", round(quantile(RMSE_test, 0.25), 2)),
            RMSE_test_3rd = sprintf("%.2f", round(quantile(RMSE_test, 0.75), 2)),
            n = n(),
            Relative_change = sprintf("%.2f", round(RMSE_test_median/RMSE_train_median - 1, 2))) %>% 
  mutate(train = paste(sprintf("%.2f", RMSE_train_median), ' (', RMSE_train_1st, ';', RMSE_train_3rd ,')', sep=''),
         test = paste(sprintf("%.2f",RMSE_test_median), ' (', RMSE_test_1st, ';', RMSE_test_3rd ,')', sep='')) %>% 
  # arrange(Relative_change)
  arrange(sample_size, num_cov, desc(Algorithm))

print(xtable(table1[, c(1,2, 13)]), include.rownames=FALSE)

################################################################################
# Mean number of terminal nodes
################################################################################

load('Results_Friedman_terminal_nodes.RData')
tab = save_results_friedman
tab$Algorithm = factor(tab$Algorithm, levels=c('CVGLMnet', 'RF', 'GB', 'BART default', 'BART', 'LM BART'), labels=c('Lasso', 'RF','GB', 'BART (200 trees)', 'BART (10 trees)', 'MOTR BART'))
tab$sample_size = factor(sprintf("%02d", as.numeric(str_split_fixed(tab$dataset, '_', n = 5)[,3])), levels=c('200', '500', '1000'), labels=c('200', '500', '1000'))
tab$num_cov = factor(sprintf("%02d", as.numeric(str_split_fixed(tab$dataset, '_', n = 5)[,5])), level=c('05', '10', '50'), labels=c('5', '10', '50'))

tab2 = tab %>%
  filter(Algorithm %in% c('MOTR BART', 'BART (200 trees)', 'BART (10 trees)')) %>% 
  group_by(Algorithm, dataset, terminal) %>% 
  mutate(num_trees = n()/10)
tab2$terminal = ifelse(tab2$terminal == -1, 'Terminal', 'Internal')

tab3 = tab2 %>% dcast(num_cov + sample_size + rep + Tree + Algorithm + dataset + num_iter + num_trees ~ terminal, value.var = "tot_par", fun.aggregate = sum)

tab4 <- tab3 %>% 
  group_by(sample_size, num_cov, rep, dataset, Algorithm, num_iter, num_trees) %>% 
  summarise(terminal_tot = round(sum(Terminal),2),
            internal_tot = round(sum(Internal),2))

tab5 = tab4 %>% 
  group_by(dataset, Algorithm, sample_size, num_cov) %>% 
  summarise(terminal_mean = formatC(mean(terminal_tot), format="d", big.mark=","),
            sd_terminal = formatC(sd(terminal_tot), format="d", big.mark=","),
            internal_mean = formatC(mean(internal_tot), format="d", big.mark=","),
            sd_internal = formatC(sd(internal_tot), format="d", big.mark=",")) %>%
  arrange(sample_size,num_cov, desc(Algorithm))

print(xtable(tab5[,c('Algorithm', 'num_cov', 'terminal_mean', 'sd_terminal')]), include.rownames=FALSE)
print(xtable(tab5[,c('Algorithm', 'num_cov', 'internal_mean', 'sd_internal')]), include.rownames=FALSE)


####################################################################################
# Real data (rank dataset) 
####################################################################################

# setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_real_DB_10_trees_1000_burn_5000_post/")
setwd("~/R/Splines_BART_Real_datasets/TEST2_Results_real_DB_10_trees_1000_burn_5000_post")
load('Results_real_DBs.RData')
# load('CLASS_Results_real_DBs.RData')

tab = save_results_real
tab$Algorithm = factor(tab$Algorithm, levels=c('CVGLMnet', 'RF', 'GB', 'BART', 'BART_def', 'LM_BART'), labels=c('Lasso', 'RF', 'GB', 'BART (10 trees)', 'BART (200 trees)', 'MOTR BART'))
tab$obs = 1:(dim(tab))[1]
head(tab)
tab$num_trees_class <- factor(tab$num_trees, levels = c('10', '200', '<NA>'), labels = c('T = 10', 'T = 200', 'No trees'))

tab_rank = tab %>% 
  group_by(Dataset, Algorithm) %>% 
  summarise(median = sprintf("%.2f", round(median(RMSE_test),2)),
            RMSE_1st_q = sprintf("%.2f", round(quantile(RMSE_test, 0.25), 2)),
            RMSE_3rd_q = sprintf("%.2f", round(quantile(RMSE_test, 0.75), 2)))%>% 
  mutate(final = paste(median, ' (', RMSE_1st_q,';' , RMSE_3rd_q, ')', sep='')) %>% 
  filter(Dataset %in% c('ankara', 'boston', 'ozone', 'compactiv', 'triazine')) %>% 
  arrange(Dataset, desc(Algorithm))

print(xtable(tab_rank[,c(1,2,6)]), include.rownames=FALSE)

####################################################################################
# Real data (number of terminal nodes) 
####################################################################################

# setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_real_DB_10_trees_1000_burn_5000_post/")
setwd("~/R/Splines_BART_Real_datasets/TEST2_Results_real_DB_10_trees_1000_burn_5000_post")
load('Results_Real_data_terminal_nodes.RData')
tab = save_real_data
tab$Algorithm = factor(tab$Algorithm, levels=c('CVGLMnet', 'RF', 'GB', 'BART default', 'BART', 'LM BART'), labels=c('Lasso', 'RF','GB', 'BART (200 trees)', 'BART (10 trees)', 'MOTR BART'))

tab2 = tab %>%
  filter(dataset %in% c('ankara', 'boston', 'ozone', 'compactiv'),
         Algorithm %in% c('MOTR BART', 'BART (200 trees)', 'BART (10 trees)')) %>% 
  group_by(Algorithm, dataset, terminal) %>% 
  mutate(num_trees = n()/10)
tab2$terminal = ifelse(tab2$terminal == -1, 'Terminal', 'Internal')

tab3 = tab2 %>% dcast(rep + Tree + Algorithm + dataset + num_iter + num_trees ~ terminal, value.var = "tot_par", fun.aggregate = sum)

tab4 <- tab3 %>% 
  group_by(rep, dataset, Algorithm, num_iter, num_trees) %>% 
  summarise(terminal_tot = round(sum(Terminal),2),
            internal_tot = round(sum(Internal),2)) %>% 
  mutate(mean_terminal_per_partition = terminal_tot/10,
         mean_internal_per_partition = internal_tot/10,
         mean_terminal_per_tree = ((terminal_tot/10)/num_iter)/num_trees,
         mean_internal_per_tree = ((internal_tot/10)/num_iter)/num_trees)

tab5 = tab4 %>% 
  group_by(dataset, Algorithm) %>% 
  summarise(terminal_mean = formatC(mean(terminal_tot), format="d", big.mark=","),
            sd_terminal = formatC(sd(terminal_tot), format="d", big.mark=","),
            internal_mean = formatC(mean(internal_tot), format="d", big.mark=","),
            sd_internal = formatC(sd(internal_tot), format="d", big.mark=",")) %>%
  arrange(dataset, desc(Algorithm))

print(xtable(tab5[,c('dataset', 'Algorithm', 'terminal_mean', 'sd_terminal')]), include.rownames=FALSE)
print(xtable(tab5[,c('dataset', 'Algorithm', 'internal_mean', 'sd_internal')]), include.rownames=FALSE)
