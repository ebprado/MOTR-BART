#' MOTR-BART for continuous outcomes
#'
#' @param x The design matrix containing the covariates that will be used in the regression
#' @param y Continuous response variable
#' @param ntrees The number of trees. Default is 10.
#' @param node_min_size Set the size of the smallest node. Default is 5.
#' @param alpha The parameter of the prior for the tree. It should be a value within (0,1). Default is 0.95.
#' @param beta The parameter of the prior for the tree. It should be a value greater than or equal to 0. Default is 2.
#' @param nu The shape of the prior for the variance.
#' @param lambda The scale of the prior for the variance.
#' @param sigma2 An inital estimate of the variance.
#' @param npost Number of MCMC iterations that should be performed after the burn-in period
#' @param nburn Number of MCMC iterations for the burn-in period
#' @param nthin The amount of thinning that should considered. If it is greater than 1, then the draws every "nthin" iterations will be returned.
#'
#' @return The following objects are returned by \code{motr_bart} function:
#'
#' \code{tau_b} = tau_b_store
#'
#' \code{log_lik} = log_lik_store
#'
#' \code{y_hat} = y_hat_store*y_sd + y_mean
#'
#' \code{full_cond}= full_cond_store

#' \code{y} = y
#'
#' \code{center_X} = center
#'
#' \code{scale_X} = scale
#'
#' \code{npost} = npost
#'
#' \code{nburn} = nburn
#'
#' \code{nthin} = nthin
#'
#' \code{store_size} = store_size
#'
#' \code{ntrees} = ntrees
#'
#' \code{y_mean} = y_mean
#'
#' \code{y_sd} = y_sd
#'
#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd'
#'
#' @examples
#'
#' # Simulate a Friedman data set
#'friedman_data = function(n, num_cov, sd_error){
#'  x = matrix(runif(n*num_cov),n,num_cov)
#'  y = 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5] + rnorm(n, sd=sd_error)
#'  return(list(y = y,
#'              x = x))
#'}
#'data = friedman_data(200, 10, 1)
#'y = data$y
#'x = data$x
#'
#'# Run MOTR-BART
#'set.seed(99)
#'fit.motr.bart = motr_bart(x, y, ntrees = 10, nburn = 100, npost = 100)
#'yhat = apply(fit.motr.bart$y_hat, 2, mean)
#'plot(y, yhat); abline(0, 1)

motr_bart = function(x,
                     y,
                     ntrees = 10,
                     node_min_size = 5,
                     alpha = 0.95,
                     beta = 2,
                     nu = 3,
                     lambda = 0.1,
                     sigma2 = 1,
                     nburn = 1000,
                     npost = 1000,
                     nthin = 1) {

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

  # Extract initial values
  log_lik = 0

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  tau_b_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  log_lik_store = rep(NA, store_size)
  full_cond_store = matrix(NA, ncol = ntrees, nrow = store_size)

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(X_orig) - 1

  # Prior of the vectors beta
  tau_b = ntrees
  V = 1/tau_b
  inv_V = tau_b

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = X)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1)

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
      tau_b_store[curr] = inv_V
      y_hat_store[curr,] = predictions
      log_lik_store[curr] = log_lik
    }

    # Start looping through trees
    for (j in 1:ntrees) {

      # Calculate partial residuals for current tree
      if(ntrees > 1) {
        current_partial_residuals = y_scale -
          get_predictions(curr_trees[-j], X, single_tree = ntrees == 2)
      } else {
        current_partial_residuals = y_scale
      }

      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune', 'change', 'swap'), 1)
      if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

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
      if((i > nburn) & (((i - nburn) %% nthin) == 0)) {
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
    predictions = get_predictions(curr_trees, X, single_tree = ntrees == 1)

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
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              store_size = store_size,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd))

} # End main function