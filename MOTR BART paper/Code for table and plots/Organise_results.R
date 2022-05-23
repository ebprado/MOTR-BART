# ---------------------------------------------------------------- #
# Description: Organise the results from simulated and real data   #
# Author: Estevao Prado                                            #
# Last modification: 14/05/2020                                    #
# -----------------------------------------------------------------#
setTimeLimit(cpu = Inf, elapsed = Inf)
library(dplyr)
library(stringr)
options(warn=0)
source('~/R/Splines_BART_Real_datasets/Spline_BART_v3.R')
source('~/R/Splines_BART_Real_datasets/Aux_functions.R')

# --------------------------------------------------------------------------------------------------------
# Real datasets
# --------------------------------------------------------------------------------------------------------

## Open the database

# setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_real_DB")
# setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_real_DB_10_trees_1000_burn_5000_post/")
# setwd("~/R/Splines_BART_Real_datasets/TEST_Results_real_DB_10_trees_1000_burn_5000_post")
setwd("~/R/Splines_BART_Real_datasets/TEST2_Results_real_DB_10_trees_1000_burn_5000_post")

db01 = LoadData('ankara');     save(db01, file='ankara_organise_results.RData')
db03 = LoadData('boston');     save(db03, file='boston_organise_results.RData')
db06 = LoadData('ozone');      save(db06, file='ozone_organise_results.RData')
db07 = LoadData('triazine');   save(db07, file='triazine_organise_results.RData')
db08 = LoadData('wine_red');   save(db08, file='wine_red_organise_results.RData')
db09 = LoadData('wine_white'); save(db09, file='wine_white_organise_results.RData')
db10 = LoadData('compactiv');  save(db10, file='compactiv_organise_results.RData')
db11 = LoadData('pole');       save(db11, file='pole_organise_results.RData')

load('ankara_organise_results.RData')
load('boston_organise_results.RData')
load('ozone_organise_results.RData')
load('triazine_organise_results.RData')
load('wine_red_organise_results.RData')
load('wine_white_organise_results.RData')
load('compactiv_organise_results.RData')
load('pole_organise_results.RData')

# save_results_real = rbind(db01,db02,db03,db04,db05,db06,db07,db08,db09)

save_results_real = rbind(db01,db03,db06,db07,db08,db09,db10,db11)

save(save_results_real, file='Results_real_DBs.RData')

# --------------------------------------------------------------------------------------------------------
# Friedman datasets
# --------------------------------------------------------------------------------------------------------

setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_Friedman_DBs/Simulation_results/")
# setwd("/users/research/eprado/R/Splines_BART_Real_datasets/Results_Friedman_DBs_10_trees_1000_burn_5000_post/Simulation_results/")

db01 = LoadData('friedman_n_200_p_5');   save(db01, file='friedman_n_200_p_5_organise_results.RData')
db02 = LoadData('friedman_n_200_p_10');  save(db02, file='friedman_n_200_p_10_organise_results.RData')
db03 = LoadData('friedman_n_200_p_50');  save(db03, file='friedman_n_200_p_50_organise_results.RData')
db04 = LoadData('friedman_n_500_p_5');   save(db04, file='friedman_n_500_p_5_organise_results.RData')
db05 = LoadData('friedman_n_500_p_10');  save(db05, file='friedman_n_500_p_10_organise_results.RData')
db06 = LoadData('friedman_n_500_p_50');  save(db06, file='friedman_n_500_p_50_organise_results.RData')
db07 = LoadData('friedman_n_1000_p_5');  save(db07, file='friedman_n_1000_p_5_organise_results.RData')
db08 = LoadData('friedman_n_1000_p_10'); save(db08, file='friedman_n_1000_p_10_organise_results.RData')
db09 = LoadData('friedman_n_1000_p_50'); save(db09, file='friedman_n_1000_p_50_organise_results.RData')

load('friedman_n_200_p_5_organise_results.RData')
load('friedman_n_200_p_10_organise_results.RData')
load('friedman_n_200_p_50_organise_results.RData')
load('friedman_n_500_p_5_organise_results.RData')
load('friedman_n_500_p_10_organise_results.RData')
load('friedman_n_500_p_50_organise_results.RData')
load('friedman_n_1000_p_5_organise_results.RData')
load('friedman_n_1000_p_10_organise_results.RData')
load('friedman_n_1000_p_50_organise_results.RData')

save_results_friedman = rbind(db01,db02,db03,db04,db05,db06,db07,db08,db09)
save(save_results_friedman, file='Results_Friedman_DBs.RData')
