# -------------------------------------------------------------------------#
# Description: it takes a BART object and predicts from it                 #
# -------------------------------------------------------------------------#

#' Title
#'
#' @param newdata The design matrix that will be used to predict the continuous outcome for.
#' @param object The fit object returned by the function \code{motr_bart}.
#' @param type Should the predicted values be the \code{"mean"} or \code{"median"} of the posterior draws? It is also possible to set \code{"all"}, where no function is applied and the posteriors draws are returned for each MCMC iteration.
#'
#' @return If \code{type = "mean"} or \code{type = "median"}, then a vector with predicted values is returned. If \code{type = "all"}, an n times npost is returned, where n is the number of row of \code{newdata}.
#' @export
#'
#' @examples
predict_motr_bart = function(object, newdata,
                         type = c('all', 'median', 'mean')) {
  # Get the means and sds to standardise the covariates from the test data
  center = object$center_X
  scale = object$scale_X
  newdata = as.matrix(cbind(1,scale(newdata, center=center, scale=scale)))

  # Create holder for predicted values
  n_newX = dim(newdata)[1]
  n_its = object$iter
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
               all = object$y_mean + object$y_sd * y_hat_mat,
               mean = object$y_mean + object$y_sd * apply(y_hat_mat,2,'mean'),
               median = object$y_mean + object$y_sd * apply(y_hat_mat,2,'median'))

  return(out)

} # end of predict function
