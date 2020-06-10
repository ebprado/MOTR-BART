#------------------------------------------------------------------------------#
# Author: Estevao Prado                                                        #
# Description code: Run MOTR BART in different real data sets.                #
# Last modified date: 09/06/2020                                               #
#------------------------------------------------------------------------------#
options(warn=2)
library(dbarts)
library(cgwtools)

save_files = '~/R/real_data_results/'
source('~/R/MOTR_BART.R')

boston     = read.csv('r_boston.csv')
ozone      = read.csv('r_ozone.csv')
triazine   = read.csv('r_triazine.csv')
ankara     = read.csv('r_ankara.csv')
wine_red   = read.csv('r_wine_red.csv')
wine_white = read.csv('r_wine_white.csv')
compactiv  = read.csv('r_compactiv.csv')
pole       = read.csv('r_pole.csv')

data_names = c("boston", "ozone", "triazine", "ankara", 'wine_red', 'wine_white', "compactiv", "pole")

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
num_trees = 10
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
    
    ### MOTR BART ---------------------------------------------------------------------------------------------
    
    num_seed = num_seed + 1
    print(list(Dataset = i,
               Repetition = j,
               seed = num_seed))
    
    set.seed(num_seed)
    
    tryCatch({
      time = system.time({
        MOTR_BART_out = rBART(X_train, y_train, num_trees = num_trees,
                            control = list(node_min_size = 5),
                            MCMC = list(iter = num_iter,
                                        burn = burn_in,
                                        thin = 1))
      })
      
    },error = function(e) e)
    
    ### BART ---------------------------------------------------------------------------------------------
    set.seed(num_seed)
    
    tryCatch({
      time = system.time({
        OrigBART = bart(X_train, y_train, X_test, ntree = num_trees, nskip = burn_in, ndpost = num_iter, keeptrees = TRUE, verbose=FALSE)
        OrigBART$my_get_trees = OrigBART$fit$getTrees()
      })
      
    }, error = function(e) e)
    
    DBtrain_file_name = paste('DBtrain', data_names[i], j, sep = '_')
    DBtest_file_name  = paste('DBtest',  data_names[i], j, sep = '_')
    LM_BART_file_name = paste('MOTR_BART', data_names[i], j, sep = '_')
    BART_file_name    = paste('BART',    data_names[i], j, sep = '_')
    
    save(db_train,      file = paste(save_files, DBtrain_file_name, '.RData', sep=''))
    save(db_test,       file = paste(save_files, DBtest_file_name,  '.RData', sep=''))
    save(MOTR_BART_out, file = paste(save_files, LM_BART_file_name, '.RData', sep=''))
    save(OrigBART,      file = paste(save_files, BART_file_name,    '.RData', sep=''))
    
  }
}
