#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'contrasts' 'model.matrix'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'

# x
# y
# sparse = TRUE
# ntrees = 10
# node_min_size = 5
# alpha = 0.95
# beta = 2
# nu = 3
# lambda = 0.1
# mu_mu = 0
# sigma2 = 1
# sigma2_mu = 1
# nburn = 1000
# npost = 1000
# nthin = 1

bart = function(   x,
                   y,
                   sparse = TRUE,
                   ntrees = 10,
                   node_min_size = 5,
                   alpha = 0.95,
                   beta = 2,
                   nu = 3,
                   lambda = 0.1,
                   mu_mu = 0,
                   sigma2 = 1,
                   sigma2_mu = 1,
                   nburn = 1000,
                   npost = 1000,
                   nthin = 1) {

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = x)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  y_hat = get_predictions(curr_trees, x, single_tree = ntrees == 1)

  # Initialise the values of the log marginalised likelihood for the stumps
  l_old = rep(tree_full_conditional(curr_trees[[1]],
                                x,
                                y_scale,
                                sigma2,
                                sigma2_mu) +
    get_tree_prior(curr_trees[[1]], alpha, beta), ntrees)


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
      y_hat_store[curr,] = y_hat
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
    }

      # Start looping through trees
      for (j in 1:ntrees) {

        current_partial_residuals = y_scale - y_hat + tree_fits_store[,j]

        # Propose a new tree via grow/change/prune/swap
        # type = sample(c('grow', 'prune', 'change', 'swap'), 1)
        type = sample(c('grow', 'prune', 'change'), 1)
        if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

        # Generate a new tree based on the current
        new_trees[[j]] = update_tree(y = y_scale,
                                     X = x,
                                     type = type,
                                     curr_tree = curr_trees[[j]],
                                     node_min_size = node_min_size,
                                     s = s)

        # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
        l_new = tree_full_conditional(new_trees[[j]],
                                      x,
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu) +
          get_tree_prior(new_trees[[j]], alpha, beta)

        # Exponentiate the results above
        a = exp(l_new - l_old[j])

        # The current tree "becomes" the new tree, if the latter is better

        if(a > runif(1)) {
          l_old[j] = l_new # the new proposed (and accepted) tree becomes the current tree

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
        curr_trees[[j]] = simulate_mu(curr_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu)
      # Updating BART predictions
      current_fit = get_predictions(curr_trees[j], x, single_tree = TRUE)
      y_hat = y_hat - tree_fits_store[,j] # subtract the old fit
      y_hat = y_hat + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

      } # End loop through trees

    sum_of_squares = sum((y_scale - y_hat)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store*y_sd^2,
              y_hat = y_hat_store*y_sd + y_mean,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              var_count_store = var_count_store,
              s = s_prob_store
              ))

} # End main function
