#--------------------------------------------------------------------------#
# Author: Estevao Prado                                                    #
# Description code: Run BART, RF, GB, GLMnet and NN for different Friedman #
#                   data sets                                              #
# Last modified date: 09/06/2020                                           #
#--------------------------------------------------------------------------#

save_files = '~/R/simulated_data_results/'

options(warn=2)
library(dbarts)
library(ranger)
library(gbm)
library(glmnet)
library(keras)
library(tensorflow)

friedman_n_200_p_5   = read.table('friedman_n_200_p_5.txt')
friedman_n_200_p_10  = read.table('friedman_n_200_p_10.txt')
friedman_n_200_p_50  = read.table('friedman_n_200_p_50.txt')
friedman_n_500_p_5   = read.table('friedman_n_500_p_5.txt')
friedman_n_500_p_10  = read.table('friedman_n_500_p_10.txt')
friedman_n_500_p_50  = read.table('friedman_n_500_p_50.txt')
friedman_n_1000_p_5  = read.table('friedman_n_1000_p_5.txt')
friedman_n_1000_p_10 = read.table('friedman_n_1000_p_10.txt')
friedman_n_1000_p_50 = read.table('friedman_n_1000_p_50.txt')

data_names = c("friedman_n_200_p_5",
               "friedman_n_200_p_10",
               "friedman_n_200_p_50",
               "friedman_n_500_p_5",
               "friedman_n_500_p_10",
               "friedman_n_500_p_50",
               'friedman_n_1000_p_5',
               'friedman_n_1000_p_10',
               "friedman_n_1000_p_50")

n_db = length(data_names)

# Putting the databases into a list

data_list = list()
for (i in 1:length(data_names)){
  db = get(data_names[[i]])
  data_list[[i]] = db
  names(data_list)[i] = data_names[i]
}

## Running BARTS

SaveResults = NULL
num_seed = 0
burn_in = 1000
n_rep = 10
num_trees = 200
num_iter = 5000

for (i in 1:length(data_names)){
  
  for (j in 1:n_rep){
    
    set.seed(j)
    n = nrow(data_list[[i]])
    index_train = sample(n, round(n*0.8))
    n_col = ncol(data_list[[i]])
    
    X_train = data_list[[i]][index_train,1:(n_col-1)]
    X_test = data_list[[i]][-index_train,1:(n_col-1)]
    y_train = data_list[[i]][index_train,n_col]
    y_test = data_list[[i]][-index_train,n_col]
    db_train = cbind(y_train, X_train)
    db_test = cbind(y_test, X_test)
    n_col_train = dim(X_train)[2]
    
    num_seed = num_seed + 1
    print(list(Dataset = i,
               Repetition = j,
               seed = num_seed))
    
    set.seed(num_seed)
    
    ### Many algorithms ---------------------------------------------------------------------------------------------
    
    tryCatch({
      time = system.time({
        OrigBART = bart(X_train, y_train, X_test, nskip = burn_in, ntree = num_trees, ndpost = num_iter, keeptrees = TRUE, verbose = FALSE)
        OrigBART$my_get_trees = OrigBART$fit$getTrees()
      })
      

      Random_Forest = ranger(y_train ~ ., data = db_train, num.trees = num_trees)
      Gradient_boosting <- gbm(formula = y_train ~ ., data =  db_train, n.trees = num_trees, distribution = 'gaussian', interaction.depth = 3)
      GLMnet <- glmnet(as.matrix(X_train), y_train, family = 'gaussian')
      CVGLMnet <- cv.glmnet(as.matrix(X_train), y_train, type.measure = 'mse')
      
    }, error = function(e) e)

    DBtrain_file_name  = paste('DBtrain',      data_names[i], j, sep = '_')
    DBtest_file_name   = paste('DBtest',       data_names[i], j, sep = '_')
    BART_file_name     = paste('BART_default', data_names[i], j, sep = '_')
    RF_file_name       = paste('RF',           data_names[i], j, sep = '_')
    GB_file_name       = paste('GB',           data_names[i], j, sep = '_')
    GLMnet_file_name   = paste('GLMnet',       data_names[i], j, sep = '_')
    CVGLMnet_file_name = paste('CVGLMnet',     data_names[i], j, sep = '_')

    save(db_train,          file = paste(save_files, DBtrain_file_name,  '.RData', sep=''))
    save(db_test,           file = paste(save_files, DBtest_file_name,   '.RData', sep=''))
    save(OrigBART,          file = paste(save_files, BART_file_name,     '.RData', sep=''))
    save(Random_Forest,     file = paste(save_files, RF_file_name,       '.RData', sep=''))
    save(Gradient_boosting, file = paste(save_files, GB_file_name,       '.RData', sep=''))
    save(GLMnet,            file = paste(save_files, GLMnet_file_name,   '.RData', sep=''))
    save(CVGLMnet,          file = paste(save_files, CVGLMnet_file_name, '.RData', sep=''))

  }
}
