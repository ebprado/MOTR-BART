# ---------------------------------------------------------------- #
# Description: Plot the results from different algorithms          #
# Author: Estevao Prado                                            #
# Last modification: 22/04/2020                                    #
# -----------------------------------------------------------------#
options(warn=1)
library(ggplot2)
library(tidyverse)
library(reshape)
library(dplyr)

setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_Friedman_DBs/Simulation_results/")
# setwd("//users/research/eprado/R/Splines_BART_Real_datasets/Results_Friedman_DBs_10_trees_1000_burn_5000_post/Simulation_results/")
load('Results_Friedman_DBs.RData')

tab = save_results_friedman
tab$Algorithm = factor(tab$Algorithm, levels=c('CVGLMnet', 'RF', 'GB', 'BART', 'BART_def', 'LM_BART'), labels=c('Lasso', 'RF','GB', 'BART (10 trees)', 'BART (200 trees)', 'MOTR BART (10 trees)'))
tab$obs = 1:(dim(tab))[1]
tab$sample_size = factor(sprintf("%02d", as.numeric(str_split_fixed(tab$Dataset, '_', n = 5)[,3])), levels=c('200', '500', '1000'), labels=c('n = 200', 'n = 500', 'n = 1000'))
tab$num_cov = factor(sprintf("%02d", as.numeric(str_split_fixed(tab$Dataset, '_', n = 5)[,5])), level=c('05', '10', '50'), labels=c('5', '10', '50'))
head(tab)
tab$num_trees_class <- factor(tab$num_trees, levels = c('10', '200', '<NA>'), labels = c('T = 10', 'T = 200', 'No trees'))

tab <- tab %>%
  group_by(Dataset) %>%
  mutate(min_train_RMSE = min(RMSE_train),
         min_test_RMSE = min(RMSE_test)) %>% 
  mutate(RRMSE_train = RMSE_train/min_train_RMSE,
         RRMSE_test = RMSE_test/min_test_RMSE)

my_plot2 <- function(varA){
  
  # pdf(paste(SaveFigures, paste0(varA,'.pdf', sep = ''), sep = ''), width = 8, height = 6)
  
  my_plot <- tab %>%
    # dplyr::filter(Dataset == varA) %>%
    # dplyr::filter(sample_size == varA) %>%
    ggplot(aes(x = num_cov, y=RMSE_test, colour=Algorithm)) + 
    geom_boxplot() +
    labs(#title = '(a)',
          # title = paste('n = ', varA, sep=''),
         # subtitle = 'Descriptive statistics',
         colour = '',
         x = 'p',
         y = 'RMSE') +
    theme_bw(base_size = 18) + 
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          legend.position = 'bottom',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    # facet_grid(. ~ num_trees_class, scales='free')
    facet_grid(.~sample_size) +
    labs(colour='')
  
  print(my_plot)
  
  # dev.off()
}

my_plot2('n = 200') # (6x6) 'Friedman_test_data.pdf'
my_plot2('n = 200') # (6x6) 'Friedman_train_data.pdf'

# --------------------------------------------------------------------------------------------------------
# Real datasets
# --------------------------------------------------------------------------------------------------------

# setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_real_DB/")
# setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_real_DB_10_trees_1000_burn_5000_post")
# setwd("~/R/Splines_BART_Real_datasets/Results_real_DB/04_all_covariates_marg_likelihood")
# setwd("~/R/Splines_BART_Real_datasets/Results_real_DB/02_Principal_components_1_2")
# setwd("~/R/Splines_BART_Real_datasets/Results_real_DB/03_Principal_components_1_sqrt_num_trees")
# setwd("~/R/Splines_BART_Real_datasets/Results_real_DB/05_standard_BART")
# setwd("~/R/Splines_BART_Real_datasets/Results_real_DB/06_MOTR_BART_intercept")
setwd('~/R/Splines_BART_Real_datasets/TEST2_Results_real_DB_10_trees_1000_burn_5000_post/')
# setwd('~/R/Splines_BART_Real_datasets/TEST_Results_real_DB_10_trees_1000_burn_5000_post/')
load('Results_real_DBs.RData')
# load('CLASS_Results_real_DBs.RData')

tab = save_results_real
tab$Algorithm = factor(tab$Algorithm, levels=c('CVGLMnet', 'RF', 'GB', 'BART', 'BART_def', 'LM_BART'), labels=c('Lasso', 'RF', 'GB', 'BART (10 trees)', 'BART (200 trees)', 'MOTR BART (10 trees)'))
tab$obs = 1:(dim(tab))[1]
head(tab)
tab$num_trees_class <- factor(tab$num_trees, levels = c('10', '200', '<NA>'), labels = c('T = 10', 'T = 200', 'No trees'))
head(tab)

tab <- tab %>%
  group_by(Dataset) %>%
  mutate(min_train_RMSE = min(RMSE_train),
         min_test_RMSE = min(RMSE_test)) %>% 
  mutate(RRMSE_train = RMSE_train/min_train_RMSE,
         RRMSE_test = RMSE_test/min_test_RMSE)
  
  ## Plots to the paper
  
  my_plot2 <- function(SaveFigure, NameFigure, varA){
    
    if (!is.null(SaveFigure) && !is.null(NameFigure)){
      pdf(paste(SaveFigure, paste0(NameFigure,'.pdf', sep = ''), sep = ''), width = 6.1, height = 12)  
    }
    
    my_plot <- tab %>%
      # dplyr::filter(Dataset %in% varA) %>%
      dplyr::filter(!(Dataset == 'compactiv' & Algorithm == 'Lasso')) %>%
      ggplot(aes(y=RMSE_test, colour=Algorithm)) + 
      geom_boxplot() +
      labs(#title = '(a)',
        # title = paste('n = ', varA, sep=''),
        # subtitle = 'Descriptive statistics',
        # x = 'Dataset',
        y = 'RMSE') +
      theme_bw(base_size = 18) + 
      theme(plot.title = element_text(size = 20, hjust = 0.5),
            legend.position = 'bottom',
            axis.text.x=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.background = element_blank()
            ) +
      # facet_grid(. ~ num_trees_class, scales='free')
      # facet_grid(Dataset ~ ., scales='free')
      facet_wrap(~ Dataset, scales='free', ncol=1) +
      labs(colour='')
    
    print(my_plot)
    
    if (!is.null(SaveFigure) && !is.null(NameFigure)){
      dev.off()
    }
  }
  
  SaveFigures = "~/R/Splines_BART_Real_datasets/"
  #NameFigure = '00_All_covariates_marginalised_likelihood'
  #NameFigure = '00_Principal_Components'
  # NameFigure = '00_Standard_BART'
  # NameFigure = '00_MOTR_BART_intercept'
  
  my_plot2(SaveFigures, 'real_data', c('ankara', 'boston', 'ozone', 'compactiv')) # better (6x12) 'real_data'
  my_plot2(SaveFigures, NULL, c('ankara', 'boston', 'ozone', 'triazine', 'wine_red', 'wine_white', 'pole', 'compactiv'
                                )) # better (6x12) 'real_data'
  my_plot2(SaveFigures, NULL, c('iono', 'haber', 'heart', 'sonar')) # better (6x12) 'real_data'
  