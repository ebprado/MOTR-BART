#' @export
predict_gam_bart = function(object, traindata, newdata, str = c('splines', 'original', 'poly'),
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

  # Splines
  var_names = names(newdata[,-1])
  newX_splines = list()
  newX_splines[[1]] = matrix(rep(1, nrow(newdata)), ncol=1)
  df = 1
  dg = 1

  aux_scale = which(scale > 0) # Removing columns where all values are equal

  # Creating the splines

  if (str == 'splines'){
    tryCatch({
      for (i in aux_scale){
        check_error = try(bs(traindata[,i], df = df, degree=dg))
        if ('try-error' %in% class(check_error)){
          X_train_splines = bs(traindata[,i], df = 1, degree=1)
          newX_splines[[i+1]] = matrix(predict(X_train_splines, newdata[,(i+1)]), ncol = 1) # 1 knot!
          names(newX_splines)[i+1] = var_names[i]
        } else {
          X_train_splines = bs(traindata[,i], df = df, degree=dg)
          newX_splines[[i+1]] = matrix(predict(X_train_splines, newdata[,(i+1)]), ncol = df) # df knots!
          names(newX_splines)[i+1] = var_names[i]
        }
      }
    },error = function(e) e)
  }

  if (str == 'original'){
    for (h in aux_scale){
      newX_splines[[h+1]] = as.matrix(newdata[,(h+1)])
    }
  }

  if (str == 'poly'){
    tryCatch({
      for (h in aux_scale){
        check_error = try(matrix(poly(traindata[,h], degree=2, raw=TRUE), nrow=nrow(traindata)))
        if ('try-error' %in% class(check_error)){
          newX_splines[[h+1]] = as.matrix(newdata[,(h+1)])
        } else {
          x_train_splines = poly(traindata[,h], degree=2, raw=TRUE)
          newX_splines[[h+1]] = predict(x_train_splines, newdata[,(h+1)])
        }
      }
    },error = function(e) e)
  }

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    newdata,
                                    newX_splines,
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
predict_gam_bart_class = function(object, newdata,
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