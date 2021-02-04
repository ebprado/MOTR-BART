#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'

motr_bart = function(x,
                     y,
                     sparse = TRUE,
                     vars_inter_slope = TRUE,
                     ntrees = 10,
                     node_min_size = 5,
                     alpha = 0.95,
                     beta = 2,
                     nu = 3,
                     lambda = 0.1,
                     sigma2 = 1,
                     nburn = 1000,
                     npost = 1000,
                     nthin = 1,
                     ancestors = FALSE) {

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
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(X_orig)
  s = rep(1/p, p)

  # Prior for the beta vector
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
  predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1, ancestors)

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

      current_partial_residuals = y_scale - predictions + tree_fits_store[,j]

      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune', 'change', 'swap'), 1)
      if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = y_scale,
                                   X = X,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s = s)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = tree_full_conditional(curr_trees[[j]],
                            X,
                            current_partial_residuals,
                            sigma2,
                            V,
                            inv_V,
                            nu,
                            lambda,
                            tau_b,
                            ancestors) +
        get_tree_prior(curr_trees[[j]], alpha, beta)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    X,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V,
                                    nu,
                                    lambda,
                                    tau_b,
                                    ancestors) +
        get_tree_prior(new_trees[[j]], alpha, beta)

      # Exponentiate the results above
      a = exp(l_new - l_old)

      if(a > runif(1)) {

        curr_trees[[j]] = new_trees[[j]] # The current tree "becomes" the new tree, if the latter is better

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
                                    X,
                                    current_partial_residuals,
                                    sigma2,
                                    inv_V,
                                    tau_b,
                                    nu,
                                    ancestors)

      current_fit = get_predictions(curr_trees[j], X, single_tree = TRUE, ancestors)
      predictions = predictions - tree_fits_store[,j] # subtract the old fit
      predictions = predictions + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Updating the predictions (y_hat)
    # predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1, ancestors)

    sum_of_squares = sum((y_scale - predictions)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    # Update sigma2_beta0 and sigma2_beta1
    if (vars_inter_slope == 'TRUE') {
      vars_betas = update_vars_intercepts_slopes(curr_trees, ntrees, sigma2)
      V = 1/c(vars_betas$var_inter, vars_betas$var_slopes)
      inv_V = 1/V
    }

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store*y_sd^2,
              y_hat = y_hat_store*y_sd + y_mean,
              center_x = center,
              scale_x = scale,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              ancestors = ancestors,
              var_count_store = var_count_store,
              s = s_prob_store,
              vars_betas = vars_betas_store))

} # End main function

########################################################################################################
# MOTR-BART for classification
########################################################################################################

#' @export

motr_bart_class = function(x,
                     y,
                     sparse = TRUE,
                     vars_inter_slope = TRUE,
                     ntrees = 10,
                     node_min_size = 5,
                     alpha = 0.95,
                     beta = 2,
                     nu = 3,
                     lambda = 0.1,
                     sigma2 = 1,
                     nburn = 1000,
                     npost = 1000,
                     nthin = 1,
                     ancestors = FALSE) {

  X_orig = x
  X = as.matrix(cbind(1,scale(x))) # standardising the covariates and adding an intercept
  y = as.integer(as.factor(y)) -1

  aux.X = apply(X, 2, unique) # Checking how many unique values each variable has
  unique.values.X = unlist(lapply(aux.X, length))
  X[, which(unique.values.X==2)] = as.matrix(X_orig[, which(unique.values.X==2)-1]) # Keeping the binary variables as they originally are

  # Quantities needed for prediction
  center = apply(X_orig, 2, mean)
  scale = apply(X_orig, 2, sd)

  center[which(unique.values.X==2)-1] = 0
  scale[which(unique.values.X==2)-1] = 1

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
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(X_orig)
  s = rep(1/p, p)

  # Prior for the beta vector
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
  predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1, ancestors)

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
    }

    # Start looping through trees
    for (j in 1:ntrees) {

      # Calculate partial residuals for current tree
      # if(ntrees > 1) {
      #   current_partial_residuals = z -
      #     get_predictions(curr_trees[-j], X, single_tree = ntrees == 2, ancestors)
      # } else {
      #   current_partial_residuals = z
      # }

      current_partial_residuals = z - predictions + tree_fits_store[,j]

      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune', 'change', 'swap'), 1)
      if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = z,
                                   X = X,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    X,
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
                                    X,
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
                                    X,
                                    current_partial_residuals,
                                    sigma2,
                                    inv_V,
                                    tau_b,
                                    nu,
                                    ancestors)

      current_fit = get_predictions(curr_trees[j], X, single_tree = TRUE, ancestors)
      predictions = predictions - tree_fits_store[,j] # subtract the old fit
      predictions = predictions + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Updating the predictions (y_hat)
    # predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1, ancestors)

    # Update z (latent variable)
    z = update_z(y, predictions)

    sum_of_squares = sum((y_scale - pnorm(predictions))^2)

    # Update sigma2_beta0 and sigma2_beta1
    if (vars_inter_slope == 'TRUE') {
      vars_betas = update_vars_intercepts_slopes(curr_trees, ntrees, sigma2)
      V = 1/c(vars_betas$var_inter, vars_betas$var_slopes)
      inv_V = 1/V
    }

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = n, nu, lambda)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
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
              ancestors = ancestors,
              var_count_store = var_count_store,
              s = s_prob_store))

} # End main function
