#' @export
predict_gam_bart = function(object, traindata, newdata,
                         type = c('all', 'median', 'mean'), test = FALSE) {
  # Get the means and sds to standardise the covariates from the test data
  remove_intercept = object$remove_intercept
  center = object$center_x
  scale = object$scale_x
  newdata_orig = newdata
  newdata = as.matrix(cbind(1,scale(newdata, center=center, scale=scale)))

  # Create holder for predicted values
  n_newX = dim(newdata)[1]
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(newdata))
  num_tress = object$num_trees

  # Splines
  var_names = names(newdata[,-1])
  newX_splines = list()
  newX_splines[[1]] = matrix(rep(1, nrow(newdata)), ncol=1)
  df = object$df
  dg = object$dg
  str = object$str
  ancestors = object$ancestors
  aux_scale = which(scale > 0) # Removing columns where all values are equal

  # Creating the splines

  if (str == 'splines'){

    tryCatch({
      for (h in aux_scale){
        check_error = try(bs(traindata[,h], df = df, degree=dg))
        if ('try-error' %in% class(check_error) || scale[h] == 1){
          newX_splines[[h+1]] = matrix(newdata_orig[,h], ncol = 1) # 1 knot!
          newdata[,(h+1)] = newX_splines[[h+1]][,1] # Get the 1st column of the splines and put it in the design matrix (that will be used to create the splitting rules)
          names(newX_splines)[h+1] = var_names[h]
        } else {
          X_train_splines = bs(traindata[,h], df = df, degree=dg)
          center = colMeans(matrix(X_train_splines, nrow=nrow(traindata)))
          sd_cov = apply(matrix(X_train_splines, nrow=nrow(traindata)),2,sd)
          if (any(is.na(sd_cov)) == TRUE) {sd_cov[which(is.na(sd_cov))] = 1}
          newX_splines[[h+1]] = matrix(scale(predict(X_train_splines, newdata_orig[,h]), center = center, scale = sd_cov), ncol = df) # 1 knot!
          # newX_splines[[h+1]] = matrix(predict(X_train_splines, newdata_orig[,h]), ncol = df) # df knots!
          newdata[,(h+1)] = newX_splines[[h+1]][,1] # Get the 1st column of the splines and put it in the design matrix (that will be used to create the splitting rules)
          names(newX_splines)[h+1] = var_names[h]
        }
      }
    },error = function(e) e)
  }

  if (str == 'original'){
    for (h in aux_scale){
      newX_splines[[h+1]] = as.matrix(newdata[,(h+1)])
    }
  }

  if (test == TRUE){
    newX_splines[[2]] = newdata_orig
    newdata[,2] = newdata_orig[,1]
  }

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    newdata,
                                    newX_splines,
                                    single_tree = length(curr_trees) == 1,
                                    ancestors = ancestors,
                                    remove_intercept)
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
predict_gam_bart_class = function(object, newdata,
                             type = c('all', 'median', 'mean')) {
  # Get the means and sds to standardise the covariates from the test data
  center = object$center_x
  scale = object$scale_x
  newdata = as.matrix(cbind(1,scale(newdata, center=center, scale=scale)))

  # Create holder for predicted values
  remove_intercept = object$remove_intercept
  n_newX = dim(newdata)[1]
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(newdata))
  num_tress = object$num_trees
  ancestors = object$ancestors

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    newdata,
                                    single_tree = length(curr_trees) == 1,
                                    ancestors = ancestors,
                                    remove_intercept = remove_intercept)
  }

  # Sort out what to return
  out = switch(type,
               all = y_hat_mat,
               mean = apply(pnorm(y_hat_mat),2,'mean'),
               median = apply(pnorm(y_hat_mat),2,'median'))

  return(out)

} # end of predict function