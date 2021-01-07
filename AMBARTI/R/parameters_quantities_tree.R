# -------------------------------------------------------------------------#
# Description: this script contains 2 functions that are used to generate  #
#              the predictions, update variance and compute the tree prior #
#              and the marginalised likelihood                             #
# -------------------------------------------------------------------------#

# 1. simulate_mu: generates the predicted values (mu's)
# 2. updata_sigma2: updates the parameters sigma2
# 3. update_z: updates the latent variables z. This is required for MOTR-BART for classification.
# 4. get_tree_prior: returns the tree log prior score
# 5. tree_full_conditional: computes the marginalised likelihood for all nodes for a given tree
# 6. get_number_distinct_cov: counts the number of distinct covariates that are used in a tree to create the splitting rules
# 7. update_linear_component:
# Compute the full conditionals -------------------------------------------------

tree_full_conditional = function(tree, R, sigma2, sigma2_mu) {

  # Function to compute log full conditional distirbution for an individual tree
  # R is a vector of partial residuals

  # Need to calculate log complete conditional, involves a sum over terminal nodes

  # First find which rows are terminal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)

  # Get node sizes for each terminal node
  nj = tree$tree_matrix[which_terminal,'node_size']

  # Get sum of residuals and sum of residuals squared within each terminal node
  sumRsq_j = aggregate(R, by = list(tree$node_indices), function(x) sum(x^2))[,2]
  S_j = aggregate(R, by = list(tree$node_indices), sum)[,2]

  # Now calculate the log posterior
  log_post = 0.5 * ( sum(log( sigma2 / (nj*sigma2_mu + sigma2))) +
              sum( (sigma2_mu* S_j^2) / (sigma2 * (nj*sigma2_mu + sigma2))))
  return(log_post)
}


# Simulate_par -------------------------------------------------------------

simulate_mu = function(tree, R, sigma2, sigma2_mu) {

  # Simulate mu values for a given tree

  # First find which rows are terminal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)

  # Get node sizes for each terminal node
  nj = tree$tree_matrix[which_terminal,'node_size']

  # Get sum of residuals in each terminal node
  sumR = aggregate(R, by = list(tree$node_indices), sum)[,2]

  # Now calculate mu values
  mu = rnorm(length(nj) ,
             mean = (sumR / sigma2) / (nj/sigma2 + sigma2_mu),
             sd = sqrt(1/(nj/sigma2 + sigma2_mu)))

  # Wipe all the old mus out for other nodes
  tree$tree_matrix[,'mu'] = NA

  # Put in just the ones that are useful
  tree$tree_matrix[which_terminal,'mu'] = mu

  return(tree)
}

# Update sigma2 -------------------------------------------------------------

update_linear_component <- function(y, yhat_bart, x, sigma2, sigma2_psi_inv){

  y_star = y - yhat_bart
  Sigma_psi = solve(t(x)%*%x / sigma2 + sigma2_psi_inv)
  Mu_psi = Sigma_psi%*%(t(x)%*%y_star/sigma2)

  beta_hat = rmvnorm(1,
                     mean = Mu_psi,
                     sigma = Sigma_psi)

  return(t(beta_hat))
}

update_g <- function(y, yhat_bart, cov_g, estimate_e, sigma2, mu_g, sigma2_g, classes_g, ng){

  sumY_e_Mu = aggregate(y - estimate_e - yhat_bart, by = list(cov_g), sum)[,2]

  sample_g = rnorm(length(ng),
        mean = (sumY_e_Mu/sigma2 + mu_g/sigma2_g)/(ng/sigma2 + 1/sigma2_g),
        sd = 1/sqrt(ng/sigma2 + 1/sigma2_g))

  new_g = rep(NA, length(y))

  for (i in 1:length(classes_g)){
    new_g[which(classes_g[i]==cov_g)] = sample_g[i]
  }
  return(list(estimate_g = new_g,
              sample_g = sample_g))
}

update_e <- function(y, yhat_bart, cov_e, estimate_g, sigma2, mu_e, sigma2_e, classes_e, ne){

  sumY_g_Mu = aggregate(y - estimate_g - yhat_bart, by = list(cov_e), sum)[,2]

  sample_e = rnorm(length(ne),
                   mean = (sumY_g_Mu/sigma2 + mu_e/sigma2_e)/(ne/sigma2 + 1/sigma2_e),
                   sd = 1/sqrt(ne/sigma2 + 1/sigma2_e))

  new_e = rep(NA, length(y))

  for (i in 1:length(classes_e)){
    new_e[which(classes_e[i]==cov_e)] = sample_e[i]
  }
  return(list(estimate_e = new_e,
              sample_e = sample_e))
}

update_sigma2 <- function(S, n, nu, lambda){
  u = 1/rgamma(1, shape = (n + nu)/2, rate = (S + nu*lambda)/2)
  return(u)
}

update_sigma2_g <- function(S_g, n_g, a_g, b_g){
  u = 1/rgamma(1, shape=(n_g/2 + a_g), rate=(S_g/2 + b_g))
  return(u)
}

update_sigma2_e <- function(S_e, n_e, a_e, b_e){
  u = 1/rgamma(1, shape=(n_e/2 + a_e), rate=(S_e/2 + b_e))
  return(u)
}

# Update the latent variable z (MOTR-BART for classification) ---------------

update_z = function(y, prediction){

  ny0 = sum(y==0)
  ny1 = sum(y==1)
  z = rep(NA, length(y))

  z[y==0] = rtruncnorm(ny0, a = -Inf, b=0,   mean = prediction[y==0], 1)
  z[y==1] = rtruncnorm(ny1, a = 0   , b=Inf, mean = prediction[y==1], 1)

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
