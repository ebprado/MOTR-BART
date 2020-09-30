# -------------------------------------------------------------------------#
# Description: this script contains 2 functions that are used to generate  #
#              the predictions, update variance and compute the tree prior #
#              and the marginalised likelihood                             #
# -------------------------------------------------------------------------#

# 1. simulate_beta: generate the 'betas' in the linear predictors
# 2. updata_sigma2: updates the parameters sigma2
# 3. update_z: updates the latent variables z. This is required for MOTR-BART for classification.
# 4. get_tree_prior: returns the tree log prior score
# 5. tree_full_conditional: computes the marginalised likelihood for all nodes for a given tree

# Compute the full conditionals -------------------------------------------------

# tree = curr_trees[[j]]
# xsplines= X_splines
# R = current_partial_residuals
tree_full_conditional = function(tree, xsplines, R, sigma2, V, inv_V, nu, lambda, tau_b, ancestors) {

  # Select the lines that correspond to terminal and internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  which_internal = which(tree$tree_matrix[,'terminal'] == 0)

  # Get the node indices for each terminal node
  curr_X_node_indices = tree$node_indices
  unique_node_indices = unique(tree$node_indices)
  log_post = NULL

  # Get the covariates that have been used as a split
  split_vars_tree <- tree$tree_matrix[which_internal, 'split_variable']
  # lm_vars <- c(1, sort(unique(as.numeric(split_vars_tree))))
  n = length(R)

  if (ancestors == FALSE) {lm_vars <- c(1, sort(unique(as.numeric(split_vars_tree))))}
  # if (ancestors == 'all covariates') {lm_vars <- 1:ncol(X)}
  if (ancestors == TRUE) {ancestors <- get_ancestors(tree)}

  # Compute the log marginalised likelihood for each terminal node
  for(i in 1:length(unique_node_indices)) {
    if (ancestors == TRUE) {
      lm_vars = c(1, ancestors[which(ancestors[,'terminal'] == unique_node_indices[i]), 'ancestor']) # Get the corresponding ancestors of the current terminal node
    }
    X_node = as.matrix(matrix(unlist(xsplines[lm_vars]), nrow=n)[curr_X_node_indices == unique_node_indices[i],]) # this is for when lm_vars = 1
    r_node = R[curr_X_node_indices == unique_node_indices[i]]
    invV = diag(ncol(X_node))*inv_V
    Lambda_node_inv = t(X_node)%*%X_node + invV
    Lambda_node = solve(t(X_node)%*%X_node + invV)
    mu_node = Lambda_node%*%((t(X_node))%*%r_node)

    log_post[i] = -0.5 * log(V) +
      0.5*log(1/det(Lambda_node_inv)) -
      (1/(2*sigma2)) * (- t(mu_node)%*%Lambda_node_inv%*%mu_node)

  }
  return(sum(log_post))
}


# Simulate_par -------------------------------------------------------------

simulate_beta = function(tree, xsplines, R, sigma2, inv_V, tau_b, nu, ancestors) {

  # First find which rows are terminal and internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  which_internal = which(tree$tree_matrix[,'terminal'] == 0)

  # Get node indices
  curr_X_node_indices = tree$node_indices
  unique_node_indices = unique(tree$node_indices)

  # Wipe all the old parameters out for other nodes
  tree$tree_matrix[,'beta_hat'] = NA

  # Get the covariates that have been used as a split
  split_vars_tree <- tree$tree_matrix[which_internal, 'split_variable']
  #lm_vars <- c(1, sort(unique((as.numeric(split_vars_tree)))))
  if (ancestors == FALSE) {lm_vars <- c(1, sort(unique(as.numeric(split_vars_tree))))}
  # if (ancestors == 'all covariates') {lm_vars <- 1:ncol(X)}
  if (ancestors == TRUE) {ancestors <- get_ancestors(tree)}
  n = length(R)

  for(i in 1:length(unique_node_indices)) {
    if (ancestors == TRUE) {
      lm_vars = c(1, ancestors[which(ancestors[,'terminal'] == unique_node_indices[i]), 'ancestor']) # Get the corresponding ancestors of the current terminal node
    }
    X_node = as.matrix(matrix(unlist(xsplines[lm_vars]), nrow=n)[curr_X_node_indices == unique_node_indices[i],])
    invV = diag(ncol(X_node))*inv_V
    r_node = R[curr_X_node_indices == unique_node_indices[i]]
    Lambda_node = solve(t(X_node)%*%X_node + invV)

    # Generate betas  -------------------------------------------------
    beta_hat = rmvnorm(1,
                 mean = Lambda_node%*%(t(X_node)%*%r_node),
                 sigma = sigma2*Lambda_node)

    # Put in the esimates
    tree$tree_matrix[unique_node_indices[i],'beta_hat'] = paste(beta_hat, collapse = ',')
  }

  return(tree)
}

# Update sigma2 -------------------------------------------------------------

update_sigma2 <- function(S, n, nu, lambda){
  u = 1/rgamma(1, shape = (n + nu)/2, rate = (S + nu*lambda)/2)
  return(u)
}

# Update the latent variable z (MOTR-BART for classification) ---------------

update_z = function(y, prediction){

  z = rnorm(length(y), prediction, sd=1)
  z[y==0] = pmin(z[y==0], 0)
  z[y==1] = pmax(z[y==1], 0)

  return(z)
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