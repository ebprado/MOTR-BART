#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'contrasts' 'model.matrix' 'as.formula'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'

# sparse = FALSE
# skip_trees = FALSE
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
# sigma2_psi = 1
# nburn = 1000
# npost = 1000
# nthin = 1
# a_g = 1
# b_g = 1
# a_e = 1
# b_e = 1

ambarti = function(x,
                   y,
                   sparse = FALSE,
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
                   sigma2_psi = 1,
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

  ng = as.numeric(tapply(cov_g, cov_g, length)) # the number of obs within each g_i
  ne = as.numeric(tapply(cov_e, cov_e, length)) # the number of obs within each e_j
  ncov = length(ng) + length(ne)

  x <- model.matrix(~ -1 + g + e, data=x,
                     contrasts.arg=list(g=contrasts(as.factor(x$g), contrasts=F),
                                        e=contrasts(as.factor(x$e), contrasts=F)))

  ###########################################
  #### Genotype
  ###########################################
  number_geno = length(ng)
  num_comb_g = floor(number_geno/2)
  formula_g = as.formula(paste('y', "~", '.^', num_comb_g))
  x_all_iter_g <- model.matrix( formula_g, data = data.frame(y = y, x[, grepl("g", colnames(x))]))
  individual_g = (1:(number_geno + 1))
  name_all_comb_g = colnames(x_all_iter_g)
  name_all_comb_g = name_all_comb_g[-c(individual_g)] # remove the individual effects

  if (number_geno%%2 == 0){ #even
    repeated_comb_g = choose(number_geno, num_comb_g)/2
    name_all_comb_g = name_all_comb_g[-c(repeated_comb_g)] # remove some equivalent columns
  }

  x_g = matrix(NA, ncol=length(name_all_comb_g), nrow=length(y))
  colnames(x_g) = name_all_comb_g

  for (k in 1:ncol(x_g)){
    name_col_g = unlist(strsplit(name_all_comb_g[k],':'))
    x_g[,k] = apply(x[,name_col_g],1,sum)
  }

  ###########################################
  #### Environment
  ###########################################

  number_env = length(ne)
  num_comb_e = floor(number_env/2)
  formula_e = as.formula(paste('y', "~", '.^', num_comb_e))
  x_all_iter_e <- model.matrix(formula_e, data = data.frame(y = y, x[, grepl("e", colnames(x))]))
  individual_e = (1:(number_env + 1))
  repeated_comb_e = choose(number_env, num_comb_e)/2
  name_all_comb_e = colnames(x_all_iter_e)
  name_all_comb_e = name_all_comb_e[-c(individual_e)] # remove individual effects

  if (length(ne)%%2 == 0){ #even
    repeated_comb_e = choose(number_env, num_comb_e)/2
    name_all_comb_e = name_all_comb_e[-c(repeated_comb_e)] # remove some equivalent columns
  }

  x_e = matrix(NA, ncol=length(name_all_comb_e), nrow=length(y))
  colnames(x_e) = name_all_comb_e

  for (k in 1:ncol(x_e)){
    name_col_e = unlist(strsplit(name_all_comb_e[k],':'))
    x_e[,k] = apply(x[,name_col_e],1,sum)
  }
  # Put x_g and x_e into a data frame and get the column indices
  x_g_e = as.data.frame(cbind(x_g, x_e))
  ind_x_g = 1:ncol(x_g)
  ind_x_e = (ncol(x_g) + 1):ncol(x_g_e)

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  bart_store = matrix(NA, ncol = length(y), nrow = store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  beta_hat_store = matrix(0, ncol = ncol(x), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)
  p_g = ncol(x_g)
  p_e = ncol(x_e)
  s_g = rep(1/p_g, p_g)
  s_e = rep(1/p_e, p_e)
  sigma2_psi = diag(p)*sigma2_psi
  sigma2_psi_inv = solve(sigma2_psi)
  yhat_linear_comp = rep(0, length(y))
  yhat_bart = rep(0, length(y))

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  yhat_bart = get_predictions(curr_trees, x, single_tree = ntrees == 1)

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
      bart_store[curr,] = yhat_bart
      y_hat_store[curr,] = y_hat
      # var_count_store[curr,] = var_count
      # s_prob_store[curr,] = s
      beta_hat_store[curr,] = beta_hat
    }

    ##### TEST #######
    ##### TEST #######
    ##### TEST #######

      beta_hat = update_linear_component(y_scale, 0, x, sigma2, sigma2_psi_inv)
      # beta_hat = update_linear_component(y_scale, yhat_bart, x, sigma2, sigma2_psi_inv)
      yhat_linear_comp = x%*%beta_hat

    ##### TEST #######
    ##### TEST #######
    ##### TEST #######

    if (skip_trees == FALSE){

      # Start looping through trees
      for (j in 1:ntrees) {

        current_partial_residuals = y_scale - (yhat_bart + yhat_linear_comp) + tree_fits_store[,j]

        # Propose a new tree via grow/change/prune/swap
        # type = sample(c('grow', 'prune', 'change', 'swap'), 1)
        type = sample(c('grow', 'prune'), 1)
        if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

        # Generate a new tree based on the current

        if (type == 'grow' || (type=='prune' && nrow(curr_trees[[j]]$tree_matrix) == 1)){

          # Below, there are two calls because we need to add an interaction of genotype and then
          # add to the same tree an interaction of environment, otherwise we run the risk of allowing
          # confunding.

          new_trees[[j]] = update_tree(y = y_scale,
                                       X = x_g_e,
                                       type = type,
                                       curr_tree = curr_trees[[j]],
                                       node_min_size = node_min_size,
                                       s = s_g,
                                       index = ind_x_g)
          var1 = new_trees[[j]]$var

          new_trees[[j]] = update_tree(y = y_scale,
                                       X = x_g_e,
                                       type = type,
                                       curr_tree = new_trees[[j]],
                                       node_min_size = node_min_size,
                                       s = s_e,
                                       index = ind_x_e)
          var2 = new_trees[[j]]$var
        } else {

        # Our context, we can't prune a terminal node because we run the risk of removing either
        # a genotype or an environment. If we remove one of them, the predicted values
        # from the tree will be confounded with the main effect associated to the enviroment/genotype
        # removed.

          new_trees[[j]] = create_stump(num_trees = 1,
                                        y = y_scale)[[1]]
        }

        # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
        l_new = tree_full_conditional(new_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu) +
          get_tree_prior(new_trees[[j]], alpha, beta)

        # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
        l_old = tree_full_conditional(curr_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu) +
          get_tree_prior(curr_trees[[j]], alpha, beta)

        # Exponentiate the results above
        a = exp(l_new - l_old)

        # The current tree "becomes" the new tree, if the latter is better

        if(a > runif(1)) {
          curr_trees[[j]] = new_trees[[j]]

          # if (type=='grow'){
          #   var_count[c(var1,var2)] = var_count[c(var1,var2)] + 1}
          #
          # if (type=='prune'){
          #   var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] - 1 } # -1 because of the intercept in X
        }

        # Update mu whether tree accepted or not
        curr_trees[[j]] = simulate_mu(curr_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu)

      current_fit = get_predictions(curr_trees[j], x_g_e, single_tree = TRUE)
      yhat_bart = yhat_bart - tree_fits_store[,j] # subtract the old fit
      yhat_bart = yhat_bart + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

      } # End loop through trees

    }

    # beta_hat = update_linear_component(y_scale, yhat_bart, x, sigma2, sigma2_psi_inv)
    #
    # # Updating predictions from the linear component
    # yhat_linear_comp = x%*%beta_hat

    # Updating the final predictions
    y_hat = yhat_linear_comp + yhat_bart

    sum_of_squares = sum((y_scale - y_hat)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    # if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
    #   s = update_s(var_count, p, 1)
    # }
  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store*y_sd^2,
              y_hat = y_hat_store*y_sd + y_mean,
              y_hat_bart = bart_store*y_sd + y_mean,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              # var_count_store = var_count_store,
              s = s_prob_store,
              beta_hat = beta_hat_store*y_sd,
              x = x,
              x_e = x_e,
              x_g = x_g,
              x_g_e = x_g_e))

} # End main function
