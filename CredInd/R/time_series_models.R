#' Fit time series models for the data of the credit indicators.
#'
#' @description This function fits the best ARIMA model and provides a forecast for data of the credit indicators.
#' @param data Should be the data of the credit indicators returned by the function \code{\link{ReadData}}.
#' @param horizon Integer that represents the number of months that will be considered to generate the forecast based on the model results.
#'
#' @return A list of objects containing the model results, correlation between the observed and fitted values, horizon, Indicator and Index. Also, a dataset containing the fitted values and the forecast results (lower and upper limits of the 95\% confidence interval) is provided.

#' @export
#' @importFrom tibble "as_data_frame"
#' @importFrom forecast "auto.arima" "forecast"
#' @importFrom stats "cor"
#' @importFrom reshape2 "melt"
#' @importFrom dplyr "select" "filter" "%>%"
# @importFrom lubridate "months"
#'
#' @seealso \code{\link{ReadData}}, \code{\link{auto.arima}}, \code{\link{forecast}} and \code{\link{cor}}.
#' @examples
#' data0 <- ReadData(Indicator = 'Interest rate', Index = 'Total')
#' model <- FitModels(data=data0, horizon=20)
#'
# library(tibble)
# library(forecast)

FitModels <-function(data, horizon){

  if (class(data) == 'CredInd' && is.numeric(horizon)){

    data_aux <- data$data
    data_aux$variable <- 'Observed'

    # Fitting the model - ARIMA(p, d, q): Autoregressive (p), d (order of the difference) and Moving average (q)
    model <- auto.arima(data_aux$value, seasonal = TRUE)

    # Getting the prediction
    prediction <- forecast(model, round(horizon))

    # Computing the correlation
    cor_obs_fit <- cor(data_aux$value, model$fitted)

    # Fixing the date format
    data_aux$Date <- as.Date(paste('01/', as.character(data_aux$Date), sep = ''), format = '%d/%m/%Y')

    # Getting the fitted values
    aux_fit <- data.frame(Date = data_aux$Date, variable='Fitted', value = model$fitted)

    # Appending fitted values to the (observed) data
    data2 <- rbind(data_aux, aux_fit)

    # Getting the confidence interval limits
    pred_low <- prediction$lower[,2] # Lower limit of the 95% confidence interval
    pred_upp <- prediction$upper[,2] # Upper limit of the 95% confidence interval

    # Setting the x-axis given the horizon for forecast
    aux_pred <- as_data_frame(seq(data_aux[dim(data_aux)[1], 'Date'], by='month', length=horizon))
    aux_pred$low_limit <- as.numeric(pred_low)
    aux_pred$upp_limit <- as.numeric(pred_upp)
    names(aux_pred) <- c('Date', 'Lower limit', 'Upper limit')

    # Preparing the database for plotting
    aux_db0 <- aux_pred %>% melt(id.vars = c('Date'), na.rm = T)
    aux_db <- rbind(aux_db0, data2)

    final_list <- list(model = model,
                       correlation = cor_obs_fit,
                       info_plot = aux_db,
                       horizon = round(horizon),
                       Indicator = data$Indicator,
                       Index = data$Index)

    class(final_list) <- 'CredInd'

    return(final_list)
  }
  else if (class(data) != 'CredInd' && is.numeric(horizon)){
        stop(paste("ERROR: The class of the data must be 'CredInd'. Generate the data by using the function ReadData.", sep=''))
  }
  else {
        stop(paste("ERROR: Either the class of the data is not 'CredInd' or the argument horizon is not numeric!", sep=''))
  }
}

# model1 <- FitModels(data = data1, horizon=20)
# model1 <- FitModels(data = data1, horizon='20')
# model1 <- FitModels(data = rnorm(100), horizon=20)
# model1 <- FitModels(data = rnorm(100), horizon='20')

