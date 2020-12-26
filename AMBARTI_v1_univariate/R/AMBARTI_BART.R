#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'contrasts' 'model.matrix'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'

# sparse = FALSE
# ntrees = 10
# node_min_size = 5
# alpha = 0.95
# beta = 2
# nu = 3
# lambda = 0.1
# mu_mu = 0
# mu_g = 0
# mu_e = 0
# sigma2 = 1
# sigma2_mu = 1
# sigma2_g = 1
# sigma2_e = 1
# nburn = 1000
# npost = 1000
# nthin = 1
# a_g = 1
# b_g = 1
# a_e = 1
# b_e = 1

ambarti = function(x,
                   y,
                   sparse = TRUE,
                   skip_trees = FALSE,
                   ntrees = 10,
                   node_min_size = 5,
                   alpha = 0.95,
                   beta = 2,
                   nu = 3,
                   lambda = 0.1,
                   mu_mu = 0,
                   mu_g = 0,
                   mu_e = 0,
                   sigma2 = 1,
                   sigma2_mu = 1,
                   sigma2_g = 1,
                   sigma2_e = 1,
                   a_g = 1,
                   b_g = 1,
                   a_e = 1,
                   b_e = 1,
                   nburn = 1000,
                   npost = 1000,
                   nthin = 1) {

  # Extract the categories for genotype and environment

  cov_g = x[,'g']
  cov_e = x[,'e']

  classes_g = sort(unique(cov_g))
  classes_e = sort(unique(cov_e))

  ng = tapply(cov_g, cov_g, length)
  ne = tapply(cov_e, cov_e, length)

  x <- model.matrix(~ -1 + g + e, data=x,
                    contrasts.arg=list(g=contrasts(as.factor(x$g), contrasts=F),
                                       e=contrasts(as.factor(x$e), contrasts=F)))

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  sigma2_g_store = rep(NA, store_size)
  sigma2_e_store = rep(NA, store_size)
  bart_store = matrix(NA, ncol = length(y), nrow = store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  g_store = matrix(0, ncol = length(ng), nrow = store_size)
  e_store = matrix(0, ncol = length(ne), nrow = store_size)

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)
  estimate_g = rep(0, n)
  estimate_e = rep(0, n)

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = x)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  bart_predictions = get_predictions(curr_trees, x, single_tree = ntrees == 1)

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
      sigma2_g_store[curr] = sigma2_g
      sigma2_e_store[curr] = sigma2_e
      bart_store[curr,] = bart_predictions
      y_hat_store[curr,] = y_hat
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      g_store[curr,] = g_
      e_store[curr,] = e_
    }

    if (skip_trees == FALSE){

      # Start looping through trees
      for (j in 1:ntrees) {

        # Calculate partial residuals for current tree
        if(ntrees > 1) {
          current_partial_residuals = y_scale -
            (estimate_g + estimate_e + get_predictions(curr_trees[-j], x, single_tree = ntrees == 2))
        } else {
          current_partial_residuals = y_scale - (estimate_g + estimate_e)
        }

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

        # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
        l_old = tree_full_conditional(curr_trees[[j]],
                                      x,
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu) +
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
        curr_trees[[j]] = simulate_mu(curr_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu)

      } # End loop through trees

    }

    # Updating the BART predictions
    bart_predictions = get_predictions(curr_trees, x, single_tree = ntrees == 1)

    # Update the estimates of genotypes and environments
    aux_g = update_g(y_scale, bart_predictions, cov_g, estimate_e, sigma2, mu_g, sigma2_g, classes_g, ng)
    aux_e = update_e(y_scale, bart_predictions, cov_e, estimate_g, sigma2, mu_e, sigma2_e, classes_e, ne)

    estimate_g <- aux_g$estimate_g
    estimate_e <- aux_e$estimate_e

    g_ = aux_g$sample_g
    e_ = aux_e$sample_e

    y_hat = estimate_g + estimate_e + bart_predictions

    sum_of_squares = sum((y_scale - y_hat)^2)
    S_g = sum((g_- mu_g)^2)
    S_e = sum((e_- mu_e)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)
    sigma2_g = update_sigma2_g(S_g, length(classes_g), a_g, b_g)
    sigma2_e = update_sigma2_e(S_e, length(classes_e), a_e, b_e)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store,
              y_hat = y_hat_store*y_sd + y_mean,
              y_hat_bart = bart_store*y_sd + y_mean,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              var_count_store = var_count_store,
              s = s_prob_store,
              sigma2_g = sigma2_g_store,
              sigma2_e = sigma2_e_store,
              g = g_store*y_sd,
              e = e_store*y_sd))

} # End main function
