#------------------------------------------------------------------------------#
# Author: Estevao Prado                                                        #
# Description code: MOTR-BART considering only covariates that are used as     #
#                   split in the linear predictor                              #
# Reference: Estevao B Prado, Rafael A Moral, Andrew Parnell (2020). Bayesian  # 
#            Additive Regression Trees with Model Trees. ArXiv preprint ArXiv: #
#            XXXXXXXX                                                          #
# Last modified date: 09/06/2020                                               #
#------------------------------------------------------------------------------#

# Comments about the code:
# 1 - If you want to understand the general structure of this code, have a look only at the main function (line 24 to 256);
# 2 - The "auxiliar" functions are presented from line 256 on in the following order:
  # 01. create_stump: initialises the trees to a stump
  # 02. update_tree: calls the corresponding function associated to the move grow, prune, change, or swap.
  # 03. grow_tree: grows a tree
  # 04. prune_tree: prunes a tree
  # 05. change_tree: changes the splitting rule that defines a pair of terminal nodes
  # 06. swap_tree: exchanges the splitting rules that define two pair of terminal nodes
  # 07. fill_tree_details: takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in each terminal node
  # 08. tree_full_conditional: computes the marginalised likelihood for all nodes for a given tree
  # 09. get_predictions: gets the predicted values from a current set of trees
  # 10. get_tree_prior: returns the tree log prior score
  # 11. simulate_mu: generate the 'betas' in the linear predictors
  # 12. updata_sigma2: updates the parameters sigma2
  # 13. get_children: it's a function that takes a node and, if the node is terminal, returns the node. If not, returns the children and calls the function again on the children
  # 14. predict_rBART: takes a BART object and predicts from it
  # 15. resample: an auxiliar function

options(warn=2)
library("devtools")
library(truncnorm) # install_github("olafmersmann/truncnorm")
library(GIGrvg)
library(MASS)
library(Matrix)
library(mvtnorm)
library(statmod)
library(MASS)

# Main function -----------------------------------------------------------

rBART = function(X, # X is the feature matrix 
                 y, # y is the target
                 num_trees = 10, # Number of trees
                 control = list(node_min_size = 5), # Size of smallest nodes
                 priors = list(alpha = 0.95, # Prior control list
                               beta = 2,
                               nu = 3,
                               lambda = 0.1), 
                 inits = list(sigma2 = 1), # sigma2 initial value
                 MCMC = list(iter = 1000, # Number of iterations
                             burn = 1000, # Size of burn in
                             thin = 1) # Amount of thinning
) { 
  
  X_orig = X
  X = as.matrix(cbind(1,scale(X))) # standardising the covariates and adding an intercept
  
  aux.X = apply(X, 2, unique) # Checking how many unique values each variable has
  unique.values.X = unlist(lapply(aux.X, length))
  X[, which(unique.values.X==2)] = as.matrix(X_orig[, which(unique.values.X==2)-1]) # Keeping the original binary variables
  
  # Quantities needed for prediction
  center = apply(X_orig, 2, mean)
  scale = apply(X_orig, 2, sd)
  
  center[which(unique.values.X==2)-1] = 0
  scale[which(unique.values.X==2)-1] = 1
  
  # Extract control parameters
  node_min_size = control$node_min_size
  
  # Extract initial values
  sigma2 = inits$sigma2
  log_lik = 0
  
  # Extract hyper-parameters
  alpha = priors$alpha # Tree shape parameter 1
  beta = priors$beta # Tree shape parameter 2
  nu = priors$nu # Sigma2 prior (parameter 1)
  lambda = priors$lambda # Sigma2 prior (parameter 2)
  
  # Extract MCMC details
  iter = MCMC$iter # Number of iterations
  burn = MCMC$burn # Size of burn in
  thin = MCMC$thin # Amount of thinning
  TotIter = burn + iter*thin # Total of iterations
  
  # Storage containers
  store_size = iter
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  tau_b_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  log_lik_store = rep(NA, store_size)
  full_cond_store = matrix(NA, ncol = num_trees, nrow = store_size)
  
  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(X_orig) - 1
  
  # Prior of the vectors beta
  tau_b = num_trees
  V = 1/tau_b
  inv_V = tau_b
  
  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = num_trees, 
                            y = y_scale,
                            X = X)
  # Initialise the new trees as current one
  new_trees = curr_trees
  
  # Initialise the predicted values to zero
  predictions = get_predictions(curr_trees, X, nu, single_tree = num_trees == 1)
  
  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60, 
                             title = 'Running rBART...')
  
  # Start the MCMC iterations loop
  for (i in 1:TotIter) {
    
    utils::setTxtProgressBar(pb, i)
    
    # If at the right place, store everything
    if((i > burn) & ((i - burn) %% thin) == 0) {
      curr = (i - burn)/thin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      tau_b_store[curr] = inv_V
      y_hat_store[curr,] = predictions
      log_lik_store[curr] = log_lik
    }
    
    # Start looping through trees
    for (j in 1:num_trees) {
      
      # Calculate partial residuals for current tree
      if(num_trees > 1) {
        current_partial_residuals = y_scale -
          get_predictions(curr_trees[-j], X, nu, single_tree = num_trees == 2)
      } else {
        current_partial_residuals = y_scale
      }
      
      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune', 'change', 'swap'), 1)
      if(i < max(floor(0.1*burn), 10)) type = 'grow' # Grow for the first few iterations 
      
      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = y_scale,
                                   X = X,
                                   type = type, 
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size)
      
      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    X,
                                    current_partial_residuals,
                                    sigma2,
                                    V,
                                    inv_V, 
                                    nu,
                                    lambda,
                                    tau_b) + 
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
                                    tau_b) + 
        get_tree_prior(curr_trees[[j]], alpha, beta)
      
      # Exponentiate the results above
      a = exp(l_new - l_old)
      
      # Save the structure of the current tree for all iterations and trees after burn-in
      if((i > burn) & (((i - burn) %% thin) == 0)) {
        full_cond_store[curr, j] = l_old
      }
      
      # The current tree "becomes" the new tree, if the latter is better
      if(a > runif(1)) {
        curr_trees[[j]] = new_trees[[j]]
      } 
      
      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_mu(curr_trees[[j]], 
                                    X,
                                    current_partial_residuals,
                                    sigma2,
                                    inv_V,
                                    tau_b, 
                                    nu)
      
    } # End loop through trees
    
    # Updating the predictions (y_hat)  
    predictions = get_predictions(curr_trees, X, nu, single_tree = num_trees == 1)
    
    S = sum((y_scale - predictions)^2)
    
    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(S, n = length(y_scale), nu, lambda)
    
    # Get the overall log likelihood
    log_lik = sum(dnorm(y_scale, mean = predictions, sd = sqrt(sigma2), log = TRUE))
    
  } # End iterations loop
  
  cat('\n') # Make sure progress bar ends on a new line
  
  return(list(trees = tree_store,
              sigma2 = sigma2_store,
              tau_b = tau_b_store,
              log_lik = log_lik_store,
              y_hat = y_hat_store*y_sd + y_mean,
              full_cond = full_cond_store,
              y = y,
              center_X = center,
              scale_X = scale,
              iter = iter,
              burn = burn,
              thin = thin,
              store_size = store_size,
              num_trees = num_trees,
              y_mean = y_mean,
              y_sd = y_sd))
  
} # End main function

# Function to create stump ------------------------------------------------

create_stump = function(num_trees,
                        y,
                        X) {
  
  # Each tree is a list of 2 elements
  # The 2 elements are the tree matrix (8 columns), and the node indices
  # The columns of the tree matrix are:
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Split variable
  # Split value
  # mu
  # Node size
  
  # Create holder for trees
  all_trees = vector('list', length = num_trees)
  # Loop through trees
  for (j in 1:num_trees) {
    # Set up each tree to have two elements in the list as described above
    all_trees[[j]] = vector('list', length = 2)
    # Give the elements names
    names(all_trees[[j]]) = c('tree_matrix', 
                              'node_indices')
    # Create the two elements: first is a matrix
    all_trees[[j]][[1]] = matrix(NA, ncol = 8, nrow = 1)
    
    # Second is the assignment to node indices
    all_trees[[j]][[2]] = rep(1, length(y))
    
    # Create column names
    colnames(all_trees[[j]][[1]]) = c('terminal', 
                                      'child_left', 
                                      'child_right',
                                      'parent',
                                      'split_variable',
                                      'split_value',
                                      'mu',
                                      'node_size')
    
    # Set values for stump 
    all_trees[[j]][[1]][1,] = c(1, NA, NA, NA, NA, NA, paste(rep(0, ncol(X)), collapse = ',') , length(y))
    
  } # End of loop through trees
  
  return(all_trees)
  
} # End of function

# Function to update trees ------------------------------------------------

update_tree = function(y, # Target variable
                       X, # Feature matrix
                       type = c('grow',   # Grow existing tree
                                'prune',  # Prune existing tree 
                                'change', # Change existing tree - change split variable and value for an internal node
                                'swap'),  # Swap existing tree - swap splitting rules for two pairs of terminal nodes
                       curr_tree,         # The current set of trees (not required if type is stump)
                       node_min_size) {   # The minimum size of a node to grow
  
  # Each tree is a list of 2 elements
  # The 2 elements are the tree matrix (8 columns), and the node indices
  # The columns of the tree matrix are:
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Split variable
  # Split value
  # mu
  # Node size
  
  # Call the appropriate function to get the new tree
  new_tree = switch(type,
                    grow = grow_tree(X, y, curr_tree, node_min_size),
                    prune = prune_tree(X, y, curr_tree),
                    change = change_tree(X, y, curr_tree, node_min_size),
                    swap = swap_tree(X, y, curr_tree, node_min_size))
  
  # Return the new tree
  return(new_tree)
  
} # End of update_tree function

# Grow_tree function ------------------------------------------------------

grow_tree = function(X, y, curr_tree, node_min_size) {
  
  # Set up holder for new tree
  new_tree = curr_tree
  
  # Get the list of terminal nodes
  terminal_nodes = which(new_tree$tree_matrix[,'terminal'] == 1)
  
  # Find terminal node sizes
  terminal_node_size = new_tree$tree_matrix[terminal_nodes,'node_size']
  
  available_values = NULL
  max_bad_trees = 10
  count_bad_trees = 0
  bad_trees = TRUE
  
  while (bad_trees ){
    
    # Set up holder for new tree
    new_tree = curr_tree
    
    # Add two extra rows to the tree in question
    new_tree$tree_matrix = rbind(new_tree$tree_matrix,
                                 c(1, NA, NA, NA, NA, NA, NA, NA), # Make sure they're both terminal
                                 c(1, NA, NA, NA, NA, NA, NA, NA))
    
    # Choose a random terminal node to split 
    node_to_split = sample(terminal_nodes, 1, 
                           prob = as.integer(as.numeric(terminal_node_size) > node_min_size)) # Choose which node to split, set prob to zero for any nodes that are too small
    
    # Choose a split variable uniformly from all columns (the first one is the intercept)
    split_variable = sample(2:ncol(X), 1)
    
    # Alternatively follow BARTMachine and choose a split value using sample on the internal values of the available
    available_values = sort(unique(X[new_tree$node_indices == node_to_split,
                                     split_variable]))
    
    if(length(available_values) == 1){
      split_value = available_values[1]
    } else if (length(available_values) == 2){
      split_value = available_values[2]
    }  else {
      # split_value = sample(available_values[-c(1,length(available_values))], 1)
      split_value = resample(available_values[-c(1,length(available_values))])
    }
    
    curr_parent = new_tree$tree_matrix[node_to_split, 'parent'] # Make sure to keep the current parent in there. Will be NA if at the root node
    new_tree$tree_matrix[node_to_split,1:6] = c(0, # Now not temrinal
                                                nrow(new_tree$tree_matrix) - 1, # child_left is penultimate row
                                                nrow(new_tree$tree_matrix),  # child_right is penultimate row
                                                curr_parent,
                                                split_variable,
                                                split_value)
    
    #  Fill in the parents of these two nodes
    new_tree$tree_matrix[nrow(new_tree$tree_matrix),'parent'] = node_to_split 
    new_tree$tree_matrix[nrow(new_tree$tree_matrix)-1,'parent'] = node_to_split 
    
    # Now call the fill function on this tree
    new_tree = fill_tree_details(new_tree, X)
    
    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[,'node_size']) <= node_min_size)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }
    
    if(count_bad_trees == max_bad_trees) return(curr_tree)
  }
  # Return new_tree
  return(new_tree)
  
} # End of grow_tree function

# Prune_tree function -----------------------------------------------------

prune_tree = function(X, y, curr_tree) {
  
  # Create placeholder for new tree
  new_tree = curr_tree
  
  if(nrow(new_tree$tree_matrix) == 1) return(new_tree) # No point in pruning a stump!
  
  # Get the list of terminal nodes
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)
  
  # Pick a random terminal node to prune
  # ONLY PICK NODES WHERE BOTH LEFT AND RIGHT CHILD ARE TERMINAL
  bad_node_to_prune = TRUE # Assume a bad node pick
  while(bad_node_to_prune) {
    
    # Choose a random terminal node
    node_to_prune = sample(terminal_nodes, 1)
    
    # Find the parent of this terminal node
    parent_pick = as.numeric(new_tree$tree_matrix[node_to_prune, 'parent'])
    
    # Get the two children of this parent
    child_left = as.numeric(new_tree$tree_matrix[parent_pick, 'child_left'])
    child_right = as.numeric(new_tree$tree_matrix[parent_pick, 'child_right'])
    
    # See whether either are terminal
    child_left_terminal = as.numeric(new_tree$tree_matrix[child_left, 'terminal'])
    child_right_terminal = as.numeric(new_tree$tree_matrix[child_right, 'terminal'])
    
    # If both are terminal then great
    if( (child_left_terminal == 1) & (child_right_terminal == 1) ) {
      bad_node_to_prune = FALSE # Have chosen a pair of terminal nodes so exist while loop
    }
    
  }# End of bad node to prune while loop
  
  # Delete these two rows from the tree matrix
  new_tree$tree_matrix = new_tree$tree_matrix[-c(child_left,child_right),,
                                              drop = FALSE]
  # Make this node terminal again with no children or split values
  new_tree$tree_matrix[parent_pick,c('terminal',
                                     'child_left',
                                     'child_right',
                                     'split_variable',
                                     'split_value')] = c(1, NA, NA, NA, NA)
  
  # If we're back to a stump no need to call fill_tree_details
  if(nrow(new_tree$tree_matrix) == 1) {
    new_tree$node_indices = rep(1, length(y))
  } else {
    # If we've removed some nodes from the middle we need to re-number all the child_left and child_right values - the parent values will still be correct
    if(node_to_prune <= nrow(new_tree$tree_matrix)) { # Only need do this if we've removed some observations from the middle of the tree matrix
      # If you're pruning any nodes which affect parent indices further down the tree then make sure to shift the parent values
      bad_parents = which(as.numeric(new_tree$tree_matrix[,'parent'])>=node_to_prune)
      # Shift them back because you have removed two rows
      new_tree$tree_matrix[bad_parents,'parent'] = as.numeric(new_tree$tree_matrix[bad_parents,'parent']) - 2
      
      for(j in node_to_prune:nrow(new_tree$tree_matrix)) {
        # Find the current parent
        curr_parent = as.numeric(new_tree$tree_matrix[j,'parent'])
        # Find both the children of this node
        curr_children = which(as.numeric(new_tree$tree_matrix[,'parent']) == curr_parent)
        # Input these children back into the parent
        new_tree$tree_matrix[curr_parent,c('child_left','child_right')] = sort(curr_children)
      } # End for loop of correcting parents and children
    } # End if statement to fill in tree details
    
    # Call the fill function on this tree
    new_tree = fill_tree_details(new_tree, X)
    
  }
  
  # Return new_tree
  return(new_tree)
  
} # End of prune_tree function

# change_tree function ----------------------------------------------------

change_tree = function(X, y, curr_tree, node_min_size) {
  
  # Change a node means change out the split value and split variable of an internal node. Need to make sure that this does now produce a bad tree (i.e. zero terminal nodes)
  
  # If current tree is a stump nothing to change
  if(nrow(curr_tree$tree_matrix) == 1) return(curr_tree)
  
  # Create a holder for the new tree
  new_tree = curr_tree
  
  # Need to get the internal nodes
  internal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 0)
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)
  
  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  max_bad_trees = 2
  count_bad_trees = 0
  bad_trees = TRUE
  while(bad_trees) {
    # Re-set the tree
    new_tree = curr_tree
    
    # choose an internal node to change
    node_to_change = sample(internal_nodes, 1)
    
    # Use the get_children function to get all the children of this node
    all_children = get_children(new_tree$tree_matrix, node_to_change)
    
    # Now find all the nodes which match these children
    use_node_indices = !is.na(match(new_tree$node_indices, all_children))
    
    # Create new split variable and value based on ignorance
    # then check this doesn't give a bad tree
    
    available_values = NULL
    
    new_split_variable = sample(2:ncol(X), 1)
    
    available_values = sort(unique(X[use_node_indices,
                                     new_split_variable]))
    
    if (length(available_values) == 1){
      new_split_value = available_values[1]
    } else if (length(available_values) == 2){
      new_split_value = available_values[2]
    } else {
      # new_split_value = sample(available_values[-c(1,length(available_values))], 1)
      new_split_value = resample(available_values[-c(1,length(available_values))])
    }
    # Update the tree details
    new_tree$tree_matrix[node_to_change,
                         c('split_variable',
                           'split_value')] = c(new_split_variable, 
                                               new_split_value)
    
    # Update the tree node indices
    new_tree = fill_tree_details(new_tree, X)
    
    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[terminal_nodes, 'node_size']) <= node_min_size)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }
    if(count_bad_trees == max_bad_trees) return(curr_tree)
    
  } # end of while loop
  
  # Return new_tree
  return(new_tree)
  
} # End of change_tree function

# swap_tree function ------------------------------------------------------

swap_tree = function(X, y, curr_tree, node_min_size) {
  
  # Swap takes two neighbouring internal nodes and swaps around their split values and variables
  
  # If current tree is a stump nothing to change
  if(nrow(curr_tree$tree_matrix) == 1) return(curr_tree)
  
  # Create a holder for the new tree
  new_tree = curr_tree
  
  # Need to get the internal nodes
  internal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 0)
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)
  
  # If less than 3 internal nodes return curr_tree
  if(length(internal_nodes) < 3) return(curr_tree)
  
  # Find pairs of neighbouring internal nodes
  parent_of_internal = as.numeric(new_tree$tree_matrix[internal_nodes,'parent'])
  pairs_of_internal = cbind(internal_nodes, parent_of_internal)[-1,]
  
  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  max_bad_trees = 2
  count_bad_trees = 0
  bad_trees = TRUE
  while(bad_trees) {
    # Re-set the tree
    new_tree = curr_tree
    
    # Pick a random pair
    nodes_to_swap = sample(1:nrow(pairs_of_internal), 1)
    
    # Get the split variables and values for this pair
    swap_1_parts = as.numeric(new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,1],
                                                   c('split_variable', 'split_value')])
    swap_2_parts = as.numeric(new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,2],
                                                   c('split_variable', 'split_value')])
    
    # Update the tree details - swap them over
    new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,1],
                         c('split_variable',
                           'split_value')] = swap_2_parts
    new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,2],
                         c('split_variable',
                           'split_value')] = swap_1_parts
    
    # Update the tree node indices
    new_tree = fill_tree_details(new_tree, X)
    
    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[terminal_nodes, 'node_size']) <= node_min_size)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }
    if(count_bad_trees == max_bad_trees) return(curr_tree)
    
  } # end of while loop
  
  # Return new_tree
  return(new_tree)
  
} # End of swap_tree function
# Fill_tree_details -------------------------------------------------------

fill_tree_details = function(curr_tree, X) {
  
  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix
  
  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix
  
  # Start with dummy node indices
  node_indices = rep(1, nrow(X))
  
  # For all but the top row, find the number of observations falling into each one
  for(i in 2:nrow(tree_matrix)) {
    
    # Get the parent
    curr_parent = as.numeric(tree_matrix[i,'parent'])
    
    # Find the split variable and value of the parent
    split_var = as.numeric(tree_matrix[curr_parent,'split_variable'])
    split_val = as.numeric(tree_matrix[curr_parent, 'split_value'])
    
    # Find whether it's a left or right terminal node
    left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                           'left', 'right')
    if(left_or_right == 'left') {
      # If left use less than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] < split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] = i
    } else {
      # If right use greater than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i
    }
  } # End of loop through table
  
  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))  
  
} # End of function

# Get complete conditions -------------------------------------------------

tree_full_conditional = function(tree, X, R, sigma2, V, inv_V, nu, lambda, tau_b) {
  
  # Select the lines that correspond to terminal and internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  which_internal = which(tree$tree_matrix[,'terminal'] == 0)
  
  # Get the node indices for each terminal node
  curr_X_node_indices = tree$node_indices
  unique_node_indices = unique(tree$node_indices)
  log_post = NULL
  
  # Get the covariates that have been used as a split
  split_vars_tree <- tree$tree_matrix[which_internal, 'split_variable']
  lm_vars <- c(1, sort(unique(as.numeric(split_vars_tree))))
  p = length(lm_vars)
  inv_V = diag(p)*inv_V
  
  # Compute the log marginalised likelihood for each terminal node
  for(i in 1:length(unique_node_indices)) {
    X_node = X[curr_X_node_indices == unique_node_indices[i], lm_vars]
    r_node = R[curr_X_node_indices == unique_node_indices[i]]
    Lambda_node_inv = t(X_node)%*%X_node + inv_V
    Lambda_node = solve(t(X_node)%*%X_node + inv_V)
    mu_node = Lambda_node%*%((t(X_node))%*%r_node)
    
    log_post[i] = -0.5 * log(V) +
      0.5*log(1/det(Lambda_node_inv)) -
      (1/(2*sigma2)) * (- t(mu_node)%*%Lambda_node_inv%*%mu_node)
    
  }
  return(sum(log_post))
}

# Get predictions ---------------------------------------------------------

get_predictions = function(trees, X, nu, single_tree = FALSE) {
  
  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]
  
  # Normally trees will be a list of lists but just in case
  if(single_tree) {
    
    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      # predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
      mu = as.numeric(unlist(strsplit(trees$tree_matrix[1, 'mu'],",")))
      predictions = X%*%mu
      
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        X_node = X[curr_X_node_indices == unique_node_indices[i],]
        mu = as.numeric(unlist(strsplit(trees$tree_matrix[unique_node_indices[i], 'mu'],",")))
        predictions[curr_X_node_indices == unique_node_indices[i]] = X_node%*%mu
        
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    
    # Do a recursive call to the function 
    partial_trees = trees
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = get_predictions(trees[[1]], X, nu, single_tree = TRUE)  + 
      get_predictions(partial_trees, X, nu,
                      single_tree = length(partial_trees) == 1)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }
  
  return(predictions)
}


# Get tree priors ---------------------------------------------------------

get_tree_prior = function(tree, alpha, beta) {
  
  # Need to work out the depth of the tree
  # First find the level of each node, then the depth is the maximum of the level
  level = rep(NA, nrow(tree$tree_matrix))
  level[1] = 0 # First row always level 0
  
  # Escpae quickly if tree is just a stump
  if(nrow(tree$tree_matrix) == 1) {
    return(log(1 - alpha)) # Tree depth is 0 
  }
  
  for(i in 2:nrow(tree$tree_matrix)) {
    # Find the current parent
    curr_parent = as.numeric(tree$tree_matrix[i,'parent'])
    # This child must have a level one greater than it's current parent
    level[i] = level[curr_parent] + 1
  }
  
  # Only compute for the internal nodes
  internal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 0)
  log_prior = 0
  for(i in 1:length(internal_nodes)) {
    log_prior = log_prior + log(alpha) - beta * log(1 + level[internal_nodes[i]]) 
  } 
  # Now add on terminal nodes
  terminal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 1)
  for(i in 1:length(terminal_nodes)) {
    log_prior = log_prior + log(1 - alpha * ((1 + level[terminal_nodes[i]])^(-beta)))
  } 
  
  
  return(log_prior)
  
}

# Simulate_par -------------------------------------------------------------

simulate_mu = function(tree, X, R, sigma2, inv_V, tau_b, nu) {
  
  # First find which rows are terminal and internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  which_internal = which(tree$tree_matrix[,'terminal'] == 0)
  
  # Get node indices
  curr_X_node_indices = tree$node_indices
  unique_node_indices = unique(tree$node_indices)
  
  # Wipe all the old parameters out for other nodes
  tree$tree_matrix[,'mu'] = NA
  
  # Get the covariates that have been used as a split
  split_vars_tree <- tree$tree_matrix[which_internal, 'split_variable']
  lm_vars <- c(1, sort(unique((as.numeric(split_vars_tree)))))
  p = length(lm_vars)
  inv_V = diag(p)*inv_V
  
  for(i in 1:length(unique_node_indices)) {
    X_node = X[curr_X_node_indices == unique_node_indices[i], lm_vars] # Only variables that have been used as split
    r_node = R[curr_X_node_indices == unique_node_indices[i]]
    Lambda_node = solve(t(X_node)%*%X_node + inv_V)
    
    # Generate betas  -------------------------------------------------
    mu = rmvnorm(1,
                 mean = Lambda_node%*%(t(X_node)%*%r_node),
                 sigma = sigma2*Lambda_node)
    
    # Put in just the ones that are useful, otherwise 0.
    aux_mu = rep(0, ncol(X))
    aux_mu[lm_vars] = mu # Only variables that have been used as split
    tree$tree_matrix[unique_node_indices[i],'mu'] = paste(aux_mu, collapse = ',')
  }
  
  return(tree)
}

# Update sigma2 -------------------------------------------------------------

update_sigma2 <- function(S, n, nu, lambda){
  u = 1/rgamma(1, shape = (n + nu)/2, rate = (S + nu*lambda)/2)
  return(u)
}

# get_children ------------------------------------------------------------

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children, 
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}

# Predict function --------------------------------------------------------

predict_rBART = function(newX, rBART_posterior, 
                         type = c('all', 'median', 'mean')) {
  # Get the means and sds to standardise the covariates from the test data
  center = rBART_posterior$center_X
  scale = rBART_posterior$scale_X
  newX = as.matrix(cbind(1,scale(newX, center=center, scale=scale)))
  
  # Create holder for predicted values
  n_newX = dim(newX)[1]
  n_its = rBART_posterior$iter
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(newX))
  num_tress = rBART_posterior$num_trees
  
  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = rBART_posterior$trees[[i]]
    
    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees, 
                                    newX,
                                    nu,
                                    single_tree = length(curr_trees) == 1)
  }
  
  # Sort out what to return
  out = switch(type,
               all = rBART_posterior$y_mean + rBART_posterior$y_sd * y_hat_mat,
               mean = rBART_posterior$y_mean + rBART_posterior$y_sd * apply(y_hat_mat,2,'mean'),
               median = rBART_posterior$y_mean + rBART_posterior$y_sd * apply(y_hat_mat,2,'median'))
  
  return(out)
  
} # end of predict function

# Sample function ----------------------------------------------------------

resample <- function(x, ...) x[sample.int(length(x), size=1), ...]
