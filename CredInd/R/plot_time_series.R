#' Visualisation of the credit indicator
#'
#' @description This function provides visualisations of the credit indicator (observed data, fitted values and lower and upper limits of the 95\% confidence interval for the forecast).
#' @param data Should be the data of the credit indicators returned by the function \code{\link{FitModels}}.
#' @export
#' @importFrom dplyr "%>%"
#' @import ggplot2
#' @seealso \code{\link{ReadData}} and \code{\link{FitModels}}.
#' @examples
#' data0 <- ReadData(Indicator = 'Interest rate', Index = 'Total')
#' model <- FitModels(data=data0, horizon=20)
#' plot_FitModels(model)

# library(ggplot2)

plot_FitModels <- function(data){

  if (class(data) == 'CredInd') {

    Date = value = variable = NULL

    db <- data$info_plot

    ggplot(db, aes(x = Date, y = value)) +
      geom_line(aes(linetype=variable, color = variable)) +
      scale_linetype_manual(values=c(Fitted = 'solid', Observed='solid', `Lower limit` = 'dotted', `Upper limit`='dotted'), limits = c('Fitted', 'Observed')) +
      scale_colour_manual(values=c(Fitted = '#F8766D', Observed='black', `Lower limit`='#F8766D', `Upper limit`='#F8766D'), limits = c('Fitted', 'Observed')) +
      labs(title = paste('Time series for', data$Indicator,'-', data$Index),
           subtitle = 'Information from the Central Bank of Brazil',
           x = 'Date',
           y = data$Indicator) +
      theme(plot.title = element_text(size=18))
  }
  else {stop("ERROR: The class of the data is not 'CredInd'. Generate the data by using the function FitModels.")}
}
# plot_FitModels(model1)
# plot_FitModels(rnorm(100)
