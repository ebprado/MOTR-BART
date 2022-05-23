library(ROCR)
library(dplyr)
library(stringr)
library(ranger)
library(gbm)
library(glmnet)
# Usefull functions

RMSE = function(yResp, yPred){
  sqrt(sum((yResp - yPred)^2)/length(yResp))
}

########################################################################################################################
# Regression
########################################################################################################################

LoadData = function(dataset){
  
  # List of all files
  all_files = list.files(path = getwd())
  all_files = as.data.frame(all_files)
  
  # Filtering the set of files associated to the dataset
  all_files %>%
    filter(str_detect(all_files, pattern = dataset)) %>% 
    arrange()
  
  Consolidate_results = NULL
  
  for (i in 1:10){
    
    # Load the RData file
    
    load(paste('BART_default_', dataset, '_', i, '.RData' , sep='')); OrigBART_default = OrigBART
    load(paste('BART_',         dataset, '_', i, '.RData' , sep=''))
    load(paste('RF_',           dataset, '_', i, '.RData' , sep=''))
    load(paste('GB_',           dataset, '_', i, '.RData' , sep=''))
    load(paste('CVGLMnet_',     dataset, '_', i, '.RData' , sep=''))
    load(paste('DBtest_',       dataset, '_', i, '.RData' , sep=''))
    load(paste('DBtrain_',      dataset, '_', i, '.RData' , sep=''))
    # load(paste('LM_BART_20_trees_',      dataset, '_', i, '.RData' , sep='')); LM_BART_out_20 = LM_BART_out
    load(paste('LM_BART_',      dataset, '_', i, '.RData' , sep=''))
    # NN <- load_model_tf(paste('NN_',           dataset, '_', i , sep=''))
    
    # Original data
    y_train = db_train[,1]
    X_train = db_train[,-1]
    y_test = db_test[,1]
    X_test = db_test[,-1]
    
    # BART_default
    
    BART_default_y_train = OrigBART_default$yhat.train.mean
    BART_default_y_test = OrigBART_default$yhat.test.mean
    
    BART_default_rmse_train = RMSE(y_train, BART_default_y_train)
    BART_default_rmse_test = RMSE(y_test, BART_default_y_test)
    aux_BART_def = c(dataset, 'BART_def', BART_default_rmse_train, BART_default_rmse_test, 200,  i)
    
    # BART
    
    BART_y_train = OrigBART$yhat.train.mean
    BART_y_test = OrigBART$yhat.test.mean
    
    BART_rmse_train = RMSE(y_train, BART_y_train)
    BART_rmse_test = RMSE(y_test, BART_y_test)
    aux_BART = c(dataset, 'BART', BART_rmse_train, BART_rmse_test, 10,  i)
    
    # LM BART
    
    LM_BART_y_train = apply(LM_BART_out$y_hat, 2, mean)
    LM_BART_y_test = predict_rBART(X_test, LM_BART_out, type='mean')
    
    LM_BART_rmse_train = RMSE(y_train, LM_BART_y_train)
    LM_BART_rmse_test = RMSE(y_test, LM_BART_y_test)
    aux_LM_BART = c(dataset, 'LM_BART', LM_BART_rmse_train, LM_BART_rmse_test, 10,  i)
    
    # LM BART (20 trees)
    
    # LM_BART_y_train_20 = apply(LM_BART_out_20$y_hat, 2, mean)
    # LM_BART_y_test_20 = predict_rBART(X_test, LM_BART_out_20, type='mean')
    # 
    # LM_BART_rmse_train_20 = RMSE(y_train, LM_BART_y_train_20)
    # LM_BART_rmse_test_20 = RMSE(y_test, LM_BART_y_test_20)
    # aux_LM_BART_20 = c(dataset, 'LM_BART_20_trees', LM_BART_rmse_train_20, LM_BART_rmse_test_20, 20,  i)
    
    # RF
    
    RF_y_train = Random_Forest$predictions
    RF_y_test = predict(Random_Forest, data = X_test)$predictions
    
    RF_rmse_train = RMSE(y_train, RF_y_train)
    RF_rmse_test = RMSE(y_test, RF_y_test)
    aux_RF = c(dataset, 'RF', RF_rmse_train, RF_rmse_test, 200,  i)
    
    # GB
    
    GB_y_train = Gradient_boosting$fit
    GB_y_test = predict(Gradient_boosting, newdata = X_test, n.trees = Gradient_boosting$n.trees)
    
    GB_rmse_train = RMSE(y_train, GB_y_train)
    GB_rmse_test = RMSE(y_test, GB_y_test)
    aux_GB = c(dataset, 'GB', GB_rmse_train, GB_rmse_test, 200,  i)
    
    # CV GLMnet
    
    CVGLMnet_y_train =  predict(CVGLMnet, newx = as.matrix(X_train), s = CVGLMnet$lambda.min)
    CVGLMnet_y_test = predict(CVGLMnet, newx = as.matrix(X_test), s = CVGLMnet$lambda.min)
    
    CVGLMnet_rmse_train = RMSE(y_train, CVGLMnet_y_train)
    CVGLMnet_rmse_test = RMSE(y_test, CVGLMnet_y_test)
    aux_CVGLMnet = c(dataset, 'CVGLMnet', CVGLMnet_rmse_train, CVGLMnet_rmse_test, NA,  i)
    
    ## NN
    
    # NN_y_train <-NN %>% predict(as.matrix(X_train))
    # NN_y_test <- NN %>% predict(as.matrix(X_test))
    # 
    # NN_rmse_train = RMSE(y_train, NN_y_train)
    # NN_rmse_test = RMSE(y_test, NN_y_test)
    # aux_NN = c(dataset, 'NN', NN_rmse_train, NN_rmse_test, NA,  i)
    
    results = as.data.frame(
      rbind(
        aux_BART_def,
        aux_BART,
        aux_LM_BART,
        aux_RF,
        aux_GB,
        aux_CVGLMnet
        # aux_LM_BART_20
        # ,aux_NN
      )
    )
    
    colnames(results) = c('Dataset', 'Algorithm', 'RMSE_train', 'RMSE_test', 'num_trees', 'repetition')
    rownames(results) = NULL
    
    results$RMSE_train = as.numeric(as.character(results$RMSE_train))
    results$RMSE_test = as.numeric(as.character(results$RMSE_test))
    
    Consolidate_results = rbind(Consolidate_results, results)
  }
  return(as.data.frame(Consolidate_results))
}


########################################################################################################################
# Regression TESTS
########################################################################################################################

LoadData2 = function(dataset){
  
  Consolidate_results = NULL
  
  for (i in 1:10){
    
    # Load the RData file
    
    
    load(paste('BART_',         dataset, '_', i, '.RData' , sep=''))
    load(paste('DBtest_',       dataset, '_', i, '.RData' , sep=''))
    load(paste('DBtrain_',      dataset, '_', i, '.RData' , sep=''))
    load(paste('LM_BART_',      dataset, '_', i, '.RData' , sep=''))
    
    # Original data
    y_train = db_train[,1]
    X_train = db_train[,-1]
    y_test = db_test[,1]
    X_test = db_test[,-1]
    
    # BART
    
    BART_y_train = OrigBART$yhat.train.mean
    BART_y_test = OrigBART$yhat.test.mean
    
    BART_rmse_train = RMSE(y_train, BART_y_train)
    BART_rmse_test = RMSE(y_test, BART_y_test)
    aux_BART = c(dataset, 'BART', BART_rmse_train, BART_rmse_test, 10,  i)
    
    # LM BART
    
    LM_BART_y_train = apply(LM_BART_out$y_hat, 2, mean)
    LM_BART_y_test = predict_rBART(X_test, LM_BART_out, type='mean')
    
    LM_BART_rmse_train = RMSE(y_train, LM_BART_y_train)
    LM_BART_rmse_test = RMSE(y_test, LM_BART_y_test)
    aux_LM_BART = c(dataset, 'LM_BART', LM_BART_rmse_train, LM_BART_rmse_test, 10,  i)
    
    results = as.data.frame(
      rbind(
        aux_BART,
        aux_LM_BART
      )
    )
    
    colnames(results) = c('Dataset', 'Algorithm', 'RMSE_train', 'RMSE_test', 'num_trees', 'repetition')
    rownames(results) = NULL
    
    results$RMSE_train = as.numeric(as.character(results$RMSE_train))
    results$RMSE_test = as.numeric(as.character(results$RMSE_test))
    
    Consolidate_results = rbind(Consolidate_results, results)
  }
  return(as.data.frame(Consolidate_results))
}



####################################################################################################
# Trees
####################################################################################################

LoadInfoTrees = function(dataset){
  
  Save = NULL
  
  for (i in 1:10){
    
    # Load the RData file
    load(paste('LM_BART_',      dataset, '_', i, '.RData' , sep=''))
    load(paste('RF_',           dataset, '_', i, '.RData' , sep=''))
    load(paste('GB_',           dataset, '_', i, '.RData' , sep=''))
    load(paste('BART_default_', dataset, '_', i, '.RData' , sep='')); OrigBART_default = OrigBART
    load(paste('BART_',         dataset, '_', i, '.RData' , sep=''))
    
    for (t in 1:200){
      
      ## Random Forests
      
      aux_RF <- treeInfo(Random_Forest, tree=t) %>%
        filter(terminal == 'TRUE') %>% 
        group_by(terminal) %>% 
        summarise(n = n())  %>% 
        mutate(rep = i,
               Tree = t,
               Algorithm = 'RF',
               dataset = dataset,
               num_iter = NA,
               tot_par = n)
      
      ## Gradient Boosting
      
      aux_GB <- pretty.gbm.tree(Gradient_boosting, i.tree = t) %>% 
        filter(SplitVar == -1) %>% 
        summarise(n = n()) %>% 
        mutate(rep = i,
               Tree = t,
               Algorithm = 'GB',
               terminal = 'TRUE',
               dataset = dataset,
               num_iter = NA,
               tot_par = n)
      
      Save = rbind(Save, aux_GB, aux_RF)
    }
    
    ## BART default
    
    BART_default_trees = OrigBART_default$my_get_trees
    BART_default_trees$var = ifelse(BART_default_trees$var == -1, -1, 0)
    num_iter_def = nrow(OrigBART_default$yhat.train)
    
    tab1 <- BART_default_trees %>%
      # filter(var == -1) %>% 
      group_by(var, tree, sample) %>% 
      summarise(terminal_nodes = n())
    
    aux_BART_default <- tab1 %>% 
      group_by(var, tree) %>%
      summarise(n = sum(terminal_nodes)) %>% 
      mutate(rep = i,
             Tree = tree,
             Algorithm = 'BART default',
             terminal = var,
             dataset = dataset,
             num_iter = num_iter_def,
             tot_par = n)
    
    Save = rbind(Save, aux_BART_default[,c('rep', 'n', 'Tree', 'Algorithm', 'terminal', 'dataset', 'num_iter', 'tot_par')])
    
    ## BART (10 trees)
    
    BART_trees = OrigBART$my_get_trees
    BART_trees$var = ifelse(BART_trees$var == -1, -1, 0)
    num_iter = nrow(OrigBART$yhat.train)
    
    tab1_BART <- BART_trees %>%
      # filter(var == -1) %>%
      group_by(var, tree, sample) %>%
      summarise(terminal_nodes = n())
    
    aux_BART <- tab1_BART %>%
      group_by(var, tree) %>%
      summarise(n = sum(terminal_nodes)) %>%
      mutate(rep = i,
             Tree = tree,
             Algorithm = 'BART',
             terminal = var,
             dataset = dataset,
             num_iter = num_iter,
             tot_par = n)
    
    Save = rbind(Save, aux_BART[,c('rep', 'n', 'Tree', 'Algorithm', 'terminal', 'dataset', 'num_iter', 'tot_par')])
    
    ## LM BART
    
    SaveBART = NULL
    SaveBART2 = NULL
    
    for (k in 1:10){
      SaveBART = rbind(SaveBART, sapply(LM_BART_out$trees, function(x) sum(as.numeric(x[[k]][[1]][,1]))))
      SaveBART2 = rbind(SaveBART2, sapply(LM_BART_out$trees, function(x) sum(!is.na(unique(as.numeric(x[[k]][[1]][,5]))))))
    }
    
    SaveBART2 = SaveBART2 + 1 # this sum is because of the intercept!
    
    # Terminal nodes
    aux_LM_BART_terminal = data.frame( rep = i,
                                       n = apply(SaveBART, 1, sum),
                                       Tree = 1:nrow(SaveBART),
                                       Algorithm = 'LM BART', 
                                       terminal = -1,
                                       dataset = dataset,
                                       num_iter = ncol(SaveBART),
                                       tot_par = apply(SaveBART2 * SaveBART, 1, sum))
    
    #Internal nodes
    aux_LM_BART_internal = data.frame(rep = i,
                                      n = apply(SaveBART-1, 1, sum),
                                      Tree = 1:nrow(SaveBART),
                                      Algorithm = 'LM BART', 
                                      terminal = 0,
                                      dataset = dataset,
                                      num_iter = ncol(SaveBART),
                                      tot_par = apply(SaveBART-1, 1, sum))
    
    Save = rbind(Save, rbind(aux_LM_BART_terminal, aux_LM_BART_internal))
    print(i)
  }
  return(Save)
}

########################################################################################################################
# Classification
########################################################################################################################

RMSE_AUC = function(yResp, yPred){
  
  RMSE = sqrt(sum((yResp - yPred)^2)/length(yResp))
  pred_ROCR = ROCR::prediction(yPred, yResp)
  auc_ROCR <- ROCR::performance(pred_ROCR, measure = "auc")
  auc_ROCR <- auc_ROCR@y.values[[1]]
  
  return(c(RMSE, auc_ROCR))
}



CLASSLoadData = function(dataset){
  
  Consolidate_results = NULL
  
  for (i in 1:10){
    
    # Load the RData file
    
    load(paste('BART_default_', dataset, '_', i, '.RData' , sep='')); OrigBART_default = OrigBART
    load(paste('BART_',         dataset, '_', i, '.RData' , sep=''))
    load(paste('RF_',           dataset, '_', i, '.RData' , sep=''))
    load(paste('GB_',           dataset, '_', i, '.RData' , sep=''))
    load(paste('CVGLMnet_',     dataset, '_', i, '.RData' , sep=''))
    load(paste('DBtest_',       dataset, '_', i, '.RData' , sep=''))
    load(paste('DBtrain_',      dataset, '_', i, '.RData' , sep=''))
    # load(paste('LM_BART_20_trees_',      dataset, '_', i, '.RData' , sep='')); LM_BART_out_20 = LM_BART_out
    load(paste('LM_BART_',      dataset, '_', i, '.RData' , sep=''))
    # NN <- load_model_tf(paste('NN_',           dataset, '_', i , sep=''))
    
    # Original data
    y_train = as.integer(as.factor(db_train[,1])) - 1
    X_train = db_train[,-1]
    y_test = as.integer(as.factor(db_test[,1])) - 1
    X_test = db_test[,-1]
    
    # BART_default
    
    BART_default_y_train = apply(pnorm(OrigBART_default$yhat.train),2,mean)
    BART_default_y_test = apply(pnorm(OrigBART_default$yhat.test),2,mean)
    
    BART_default_rmse_train = RMSE_AUC(y_train, BART_default_y_train)
    BART_default_rmse_test = RMSE_AUC(y_test, BART_default_y_test)
    aux_BART_def = c(dataset, 'BART_def', BART_default_rmse_train, BART_default_rmse_test, 200,  i)
    
    # BART
    
    BART_y_train = apply(pnorm(OrigBART$yhat.train),2,mean)
    BART_y_test = apply(pnorm(OrigBART$yhat.test),2,mean)
    
    BART_rmse_train = RMSE_AUC(y_train, BART_y_train)
    BART_rmse_test = RMSE_AUC(y_test, BART_y_test)
    aux_BART = c(dataset, 'BART', BART_rmse_train, BART_rmse_test, 10,  i)
    
    # LM BART
    
    LM_BART_y_train = apply(LM_BART_out$y_hat, 2, mean)
    LM_BART_y_test = predict_rBART(X_test, LM_BART_out, type='mean')
    
    LM_BART_rmse_train = RMSE_AUC(y_train, LM_BART_y_train)
    LM_BART_rmse_test = RMSE_AUC(y_test, LM_BART_y_test)
    aux_LM_BART = c(dataset, 'LM_BART', LM_BART_rmse_train, LM_BART_rmse_test, 10,  i)
    
    # LM BART (20 trees)
    
    # LM_BART_y_train_20 = apply(LM_BART_out_20$y_hat, 2, mean)
    # LM_BART_y_test_20 = predict_rBART(X_test, LM_BART_out_20, type='mean')
    # 
    # LM_BART_rmse_train_20 = RMSE_AUC(y_train, LM_BART_y_train_20)
    # LM_BART_rmse_test_20 = RMSE_AUC(y_test, LM_BART_y_test_20)
    # aux_LM_BART_20 = c(dataset, 'LM_BART_20_trees', LM_BART_rmse_train_20, LM_BART_rmse_test_20, 20,  i)
    
    # RF
    
    RF_y_train = as.integer(as.factor(Random_Forest$predictions)) -1
    RF_y_test = as.integer(as.factor(predict(Random_Forest, data = X_test)$predictions))-1
    
    RF_rmse_train = RMSE_AUC(y_train, RF_y_train)
    RF_rmse_test = RMSE_AUC(y_test, RF_y_test)
    aux_RF = c(dataset, 'RF', RF_rmse_train, RF_rmse_test, 200,  i)
    
    # GB
    
    GB_y_train = Gradient_boosting$fit
    GB_y_test = predict(Gradient_boosting, newdata = X_test, n.trees = Gradient_boosting$n.trees)
    
    GB_rmse_train = RMSE_AUC(y_train, GB_y_train)
    GB_rmse_test = RMSE_AUC(y_test, GB_y_test)
    aux_GB = c(dataset, 'GB', GB_rmse_train, GB_rmse_test, 200,  i)
    
    # CV GLMnet
    
    CVGLMnet_y_train =  predict(CVGLMnet, newx = as.matrix(X_train), s = CVGLMnet$lambda.min)
    CVGLMnet_y_test = predict(CVGLMnet, newx = as.matrix(X_test), s = CVGLMnet$lambda.min)
    
    CVGLMnet_rmse_train = RMSE_AUC(y_train, CVGLMnet_y_train)
    CVGLMnet_rmse_test = RMSE_AUC(y_test, CVGLMnet_y_test)
    aux_CVGLMnet = c(dataset, 'CVGLMnet', CVGLMnet_rmse_train, CVGLMnet_rmse_test, NA,  i)
    
    ## NN
    
    # NN_y_train <-NN %>% predict(as.matrix(X_train))
    # NN_y_test <- NN %>% predict(as.matrix(X_test))
    # 
    # NN_rmse_train = RMSE_AUC(y_train, NN_y_train)
    # NN_rmse_test = RMSE_AUC(y_test, NN_y_test)
    # aux_NN = c(dataset, 'NN', NN_rmse_train, NN_rmse_test, NA,  i)
    
    results = as.data.frame(
      rbind(
        aux_BART_def,
        aux_BART,
        aux_LM_BART,
        aux_RF,
        aux_GB,
        aux_CVGLMnet
        # aux_LM_BART_20
        # ,aux_NN
      )
    )
    
    colnames(results) = c('Dataset', 'Algorithm', 'RMSE_train', 'AUC_train', 'RMSE_test', 'AUC_test', 'num_trees', 'repetition')
    rownames(results) = NULL
    
    results$RMSE_train = as.numeric(as.character(results$RMSE_train))
    results$RMSE_test = as.numeric(as.character(results$RMSE_test))
    results$AUC_train = as.numeric(as.character(results$AUC_train))
    results$AUC_test = as.numeric(as.character(results$AUC_test))
    
    Consolidate_results = rbind(Consolidate_results, results)
  }
  return(as.data.frame(Consolidate_results))
}

