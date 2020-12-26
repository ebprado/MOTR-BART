#' @export
predict_ambarti = function(object, newdata,
                         type = c('all', 'median', 'mean')) {
  # Get the means and sds to standardise the covariates from the test data

  center = object$center_x
  scale = object$scale_x
  newdata = as.matrix(cbind(1,scale(newdata, center=center, scale=scale)))

  # Create holder for predicted values
  n_newX = dim(newdata)[1]
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(newdata))
  num_tress = object$num_trees

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    newdata,
                                    single_tree = length(curr_trees) == 2)
  }

  # Sort out what to return
  out = switch(type,
               all = object$y_mean + object$y_sd * y_hat_mat,
               mean = object$y_mean + object$y_sd * apply(y_hat_mat,2,'mean'),
               median = object$y_mean + object$y_sd * apply(y_hat_mat,2,'median'))

  return(out)

} # end of predict function


########################################################################################################
# Predictions for classification
########################################################################################################

#' @export
predict_ambarti_class = function(object,
                                   newdata,
                                   type = c('all', 'median', 'mean')) {
  # Get the means and sds to standardise the covariates from the test data
  center = object$center_x
  scale = object$scale_x
  newdata = as.matrix(cbind(1,scale(newdata, center=center, scale=scale)))

  # Create holder for predicted values
  n_newX = dim(newdata)[1]
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(newdata))
  num_tress = object$num_trees

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    newdata,
                                    single_tree = length(curr_trees) == 1)
  }

  # Sort out what to return
  out = switch(type,
               all = y_hat_mat,
               mean = apply(pnorm(y_hat_mat),2,'mean'),
               median = apply(pnorm(y_hat_mat),2,'median'))

  return(out)

} # end of predict function
