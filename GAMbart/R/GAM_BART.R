#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'poly' 'predict'
#' @importFrom splines 'bs'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom Matrix 'bdiag'

x
y
sparse = TRUE
vars_inter_slope = TRUE
str = c('splines')
df = 10
dg = 3
ntrees = 10
node_min_size = 10
alpha = 0.95
beta = 2
nu = 3
lambda = 0.1
sigma2 = 1
nburn = 1000
npost = 1000
nthin = 1
ancestors = FALSE
one_var_per_tree = FALSE
remove_intercept = FALSE
test = FALSE
penalty = 'EM'

gam_bart = function(x,
                    y,
                    sparse = TRUE,
                    vars_inter_slope = TRUE,
                    str = c('splines', 'original'),
                    df = 10,
                    dg = 3,
                    ntrees = 10,
                    node_min_size = 10,
                    alpha = 0.95,
                    beta = 2,
                    nu = 3,
                    lambda = 0.1,
                    sigma2 = 1,
                    nburn = 1000,
                    npost = 1000,
                    nthin = 1,
                    ancestors = FALSE,
                    one_var_per_tree = FALSE,
                    remove_intercept = FALSE,
                    test = FALSE,
                    penalty = 'ridge') {

  X_orig = x
  X = as.matrix(cbind(1,x)) # adding an intercept

  aux.X = apply(X, 2, unique) # Checking how many unique values each variable has
  unique.values.X = unlist(lapply(aux.X, length))
  X[, which(unique.values.X==2)] = as.matrix(X_orig[, which(unique.values.X==2)-1]) # Keeping the binary variables as they originally are

  # Quantities needed for prediction
  center = apply(X_orig, 2, mean)
  scale = apply(X_orig, 2, sd)

  center[which(unique.values.X==2)-1] = 0
  scale[which(unique.values.X==2)-1] = 1

  aux_scale = which(scale > 0) # Removing columns where all values are equal

  var_names = names(X_orig)
  X_splines = list()
  X_splines[[1]] = matrix(rep(1, nrow(X_orig)), ncol=1)

  # Create the splines ----------------------------------------------------------------------
  if (str == 'splines'){

    tryCatch({

      for (h in aux_scale){
        check_error = try(bs(X_orig[,h], df = df, degree=dg))
        if ('try-error' %in% class(check_error) || h %in% (which(unique.values.X==2)-1)){ # binary variables
          X_splines[[h+1]] = matrix(X_orig[,h], ncol = 1) # 1 knot!
          # X[,(h+1)] = X_splines[[h+1]][,1] # Get the 1st column of the splines and put it in the design matrix (that will be used to create the splitting rules)
          names(X_splines)[h+1] = var_names[h]
        } else {
          X_splines[[h+1]] = matrix(scale(bs(X_orig[,h], df = df, degree = dg)), ncol = df) # df knots!
          sd_scaled_splines = apply(X_splines[[h+1]], 2, sd)
          if (any(is.na(sd_scaled_splines)) == TRUE) {
            X_splines[[h+1]][,which(is.na(sd_scaled_splines))] = 0
          }
          # X[,(h+1)] = X_splines[[h+1]][,1]
          names(X_splines)[h+1] = var_names[h]
        }
      }
    },error = function(e) e)
  }

  # Keep the (standardised) original covariates ------------------------------------------------------------
  if (str == 'original'){
    for (h in aux_scale){
      X_scaled = scale(X_orig) # standardising the covariates
      X_splines[[h+1]] = as.matrix(scale(X_scaled[,(h+1)]))
    }
  }

  if (test == TRUE){
    X_splines[[2]] = X_orig
    X[,2] = X_orig[,1]
  }

  # Get the number of columns/basis generated for each variable
  num_columns = rep(NA, length(X_splines) - 1) # Discount the intercept
  for (h in 1:length(num_columns)){
    num_columns[h] = ncol(X_splines[[h+1]])
  }

  # Create the penalty matrices considering the number of basis functions generated for each variable
  if (penalty == 'EM'){
    penalty_matrix = list()
    penalty_matrix[[1]] = matrix(1, ncol=1)
    for (h in 1:length(num_columns)){
      ncolumns = num_columns[h]
      P <- diff(diag(ncolumns), differences = 1)
      penalty_matrix[[h+1]] <- t(P)%*%P
    }
  }
  if (penalty == 'ridge'){
    penalty_matrix = list()
    penalty_matrix[[1]] = matrix(1, ncol=1)
    for (h in 1:length(num_columns)){
      ncolumns = num_columns[h]
      penalty_matrix[[h+1]] <- diag(ncolumns)
    }
  }

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(X_orig))
  var_count_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  vars_betas_store = matrix(0, ncol = 2, nrow = store_size)

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(X_orig)
  s = rep(1/p, p)

  # If the one_var_per_tree = TRUE, then the number of trees is the number of covariates in the data set

  if (one_var_per_tree == TRUE){
    ntrees = ncol(X_orig)
  }

  # Prior of the vectors beta
  tau_b = ntrees
  V = rep(1/tau_b, 2)
  inv_V = 1/V

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = X)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  predictions = get_predictions(curr_trees, X, X_splines, single_tree = ntrees == 1, ancestors, remove_intercept)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = predictions
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      vars_betas_store[curr,] = V
    }

    # Start looping through trees
    for (j in 1:ntrees) {

      # Calculate partial residuals for current tree
      if(ntrees > 1) {
        current_partial_residuals = y_scale -
          get_predictions(curr_trees[-j], X, X_splines, single_tree = ntrees == 2, ancestors, remove_intercept)
      } else {
        current_partial_residuals = y_scale
      }

      # Propose a new tree via grow/change/prune/swap

      if (one_var_per_tree == TRUE)
        {type = sample(c('grow', 'prune'), 1)}
      else
        {type = sample(c('grow', 'prune', 'change', 'swap'), 1)}

      if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

      # Generate a new tree based on the current

      new_trees[[j]] = update_tree(y = y_scale,
                                   X = X,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s = s,
                                   index_tree = j,
                                   one_var_per_tree = one_var_per_tree)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    X_splines,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    ancestors,
                                    remove_intercept,
                                    penalty_matrix) +
        get_tree_prior(new_trees[[j]], alpha, beta)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = tree_full_conditional(curr_trees[[j]],
                                    X_splines,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    ancestors,
                                    remove_intercept,
                                    penalty_matrix) +
        get_tree_prior(curr_trees[[j]], alpha, beta)

      # Exponentiate the results above
      a = exp(l_new - l_old)

      # The current tree "becomes" the new tree, if the latter is better
      if(a > runif(1)) {
        curr_trees[[j]] = new_trees[[j]]

        if (type =='change'){
          var_count[curr_trees[[j]]$var[1] - 1] = var_count[curr_trees[[j]]$var[1] - 1] - 1
          var_count[curr_trees[[j]]$var[2] - 1] = var_count[curr_trees[[j]]$var[2] - 1] + 1
        }

        if (type=='grow'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] + 1 } # -1 because of the intercept in X

        if (type=='prune'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] - 1 } # -1 because of the intercept in X
      }

      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_beta(curr_trees[[j]],
                                    X_splines,
                                    current_partial_residuals,
                                    sigma2,
                                    inv_V,
                                    tau_b,
                                    nu,
                                    ancestors,
                                    remove_intercept,
                                    penalty_matrix)

    } # End loop through trees

    # Updating the predictions (y_hat)
    predictions = get_predictions(curr_trees, X, X_splines, single_tree = ntrees == 1, ancestors, remove_intercept)

    S = sum((y_scale - predictions)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(S, n = length(y_scale), nu, lambda)

    # Update sigma2_beta0 and sigma2_beta1
    if (vars_inter_slope == 'TRUE') {
      vars_betas = update_vars_intercepts_slopes(curr_trees, ntrees, sigma2)
      V = 1/c(vars_betas$var_inter, vars_betas$var_slopes)
      inv_V = 1/V
    }

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE'){
      s = update_s(var_count, p, 1)
    }
  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store,
              y_hat = y_hat_store*y_sd + y_mean,
              center_x = center,
              scale_x = scale,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              df = df,
              dg = dg,
              str = str,
              ancestors = ancestors,
              var_count_store = var_count_store,
              s = s_prob_store,
              vars_betas = vars_betas_store,
              remove_intercept = remove_intercept))

} # End main function

########################################################################################################
# MOTR-BART for classification
########################################################################################################

#' @export

gam_bart_class = function(x,
                          y,
                          sparse = TRUE,
                          vars_inter_slope = TRUE,
                          str = c('splines', 'original'),
                          df = 1,
                          dg = 1,
                          ntrees = 10,
                          node_min_size = 10,
                          alpha = 0.95,
                          beta = 2,
                          nu = 3,
                          lambda = 0.1,
                          sigma2 = 1,
                          nburn = 1000,
                          npost = 1000,
                          nthin = 1,
                          ancestors = FALSE,
                          remove_intercept = FALSE) {

  X_orig = x
  X = as.matrix(cbind(1,scale(x))) # standardising the covariates and adding an intercept

  aux.X = apply(X, 2, unique) # Checking how many unique values each variable has
  unique.values.X = unlist(lapply(aux.X, length))
  X[, which(unique.values.X==2)] = as.matrix(X_orig[, which(unique.values.X==2)-1]) # Keeping the binary variables as they originally are

  # Quantities needed for prediction
  center = apply(X_orig, 2, mean)
  scale = apply(X_orig, 2, sd)

  center[which(unique.values.X==2)-1] = 0
  scale[which(unique.values.X==2)-1] = 1

  aux_scale = which(scale > 0) # Removing columns where all values are equal

  var_names = names(X_orig)
  X_splines = list()
  X_splines[[1]] = matrix(rep(1, nrow(X_orig)), ncol=1)

  # Create the splines ----------------------------------------------------------------------
  if (str == 'splines'){

    tryCatch({

      for (h in aux_scale){
        check_error = try(bs(X_orig[,h], df = df, degree=dg))
        if ('try-error' %in% class(check_error) || h %in% (which(unique.values.X==2)-1)){ # binary variables
          X_splines[[h+1]] = matrix(X_orig[,h], ncol = 1) # 1 knot!
          X[,(h+1)] = X_splines[[h+1]][,1] # Get the 1st column of the splines and put it in the design matrix (that will be used to create the splitting rules)
          names(X_splines)[h+1] = var_names[h]
        } else {
          X_splines[[h+1]] = matrix(scale(bs(X_orig[,h], df = df, degree = dg)), ncol = df) # df knots!
          X[,(h+1)] = X_splines[[h+1]][,1]
          names(X_splines)[h+1] = var_names[h]
        }
      }
    },error = function(e) e)
  }

  # Keep the (standardised) original covariates ------------------------------------------------------------
  if (str == 'original'){
    for (h in aux_scale){
      X_splines[[h+1]] = as.matrix(X[,(h+1)])
    }
  }

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(X_orig))
  var_count_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(X_orig), nrow = store_size)
  vars_betas_store = matrix(0, ncol = 2, nrow = store_size)

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(X_orig)
  s = rep(1/p, p)

  # Prior of the vectors beta
  tau_b = ntrees
  V = rep(1/tau_b, 2)
  inv_V = 1/V

  # Initial values
  z = ifelse(y == 0, -3, 3)

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = z,
                            X = X)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  predictions = get_predictions(curr_trees, X, X_splines, single_tree = ntrees == 1, ancestors, remove_intercept)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = pnorm(predictions)
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      vars_betas_store[curr,] = V
    }

    # Start looping through trees
    for (j in 1:ntrees) {

      # Calculate partial residuals for current tree
      if(ntrees > 1) {
        current_partial_residuals = z -
          get_predictions(curr_trees[-j], X, X_splines, single_tree = ntrees == 2, ancestors, remove_intercept)
      } else {
        current_partial_residuals = z
      }

      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune', 'change', 'swap'), 1)
      if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = z,
                                   X = X,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s = s)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    X_splines,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    ancestors) +
        get_tree_prior(new_trees[[j]], alpha, beta)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = tree_full_conditional(curr_trees[[j]],
                                    X_splines,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    ancestors) +
        get_tree_prior(curr_trees[[j]], alpha, beta)

      # Exponentiate the results above
      a = exp(l_new - l_old)

      # The current tree "becomes" the new tree, if the latter is better
      if(a > runif(1)) {
        curr_trees[[j]] = new_trees[[j]]

        if (type == 'change'){
          var_count[curr_trees[[j]]$var[1] - 1] = var_count[curr_trees[[j]]$var[1] - 1] - 1
          var_count[curr_trees[[j]]$var[2] - 1] = var_count[curr_trees[[j]]$var[2] - 1] + 1
        }

        if (type == 'grow'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] + 1 } # -1 because of the intercept in X

        if (type == 'prune'){
          var_count[curr_trees[[j]]$var - 1] = var_count[curr_trees[[j]]$var - 1] - 1 } # -1 because of the intercept in X
      }

      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_beta(curr_trees[[j]],
                                      X_splines,
                                      current_partial_residuals,
                                      sigma2,
                                      inv_V,
                                      tau_b,
                                      nu,
                                      ancestors)

    } # End loop through trees

    # Updating the predictions (y_hat)
    predictions = get_predictions(curr_trees, X, X_splines, single_tree = ntrees == 1, ancestors, remove_intercept)

    # Update z (latent variable)
    z = update_z(y, predictions)

    S = sum((y_scale - predictions)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(S, n = length(y_scale), nu, lambda)

    # Update sigma2_beta0 and sigma2_beta1
    if (vars_inter_slope == 'TRUE') {
      vars_betas = update_vars_intercepts_slopes(curr_trees, ntrees, sigma2)
      V = c(vars_betas$var_inter, vars_betas$var_slopes)
      inv_V = 1/V
    }

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE'){
      s = update_s(var_count, p, 1)
    }

  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store,
              y_hat = y_hat_store,
              center_x = center,
              scale_x = scale,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              str = str,
              ancestors = ancestors,
              var_count_store = var_count_store,
              s = s_prob_store,
              vars_betas = vars_betas_store,
              remove_intercept = remove_intercept))

} # End main function

