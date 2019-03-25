#' Databases of the credit indicators from the Central Bank of Brazil.
#'
#' @description This function reads and manipulates the databases of the credit indicators of the public and private banks from Brazil. The information was obtained in the site of the Central Bank of Brazil.
#' @param Indicator The credit indicator can be either: Allowance, Due loans, Due loans by capital, Spread, House debt and Interest rate.
#' @param Index Given the Indicator, the parameter Index identifies which information should be consider. The parameter Index varies depending on the Indicator (see below Details).
#' @details For instance, when Indicator = 'Allowance', the parameter Index can be Total, Public banks, National private control, Foreign control and Private control.
#'\itemize{
#'  \item{Total: }{Percentage of total provision in relation to the loan portfolio of the financial system.}
#'  \item{Public banks: }{Percentage of total provision in relation to the loan portfolio of the financial institutions under public control.}
#'  \item{National private control: }{Percentage of total provision in relation to the loan portfolio of the financial institutions under national private control.}
#'  \item{Foreign control: }{Percentage of total provision in relation to the loan portfolio of the financial institutions under foreign control.}
#'  \item{Private control: }{Percentage of total provision in relation to the loan portfolio of the financial institutions under private control.}
#'}
#'
#' @details For Indicator = 'Spread', the parameter Index can be Total, New operations companies, New operations households, Non-revolving total, Non-revolving companies and Non-revolving households.
#'\itemize{
#'  \item{Total: }{Average spread of new credit operations - Total - (percentage points) p.p.}
#'  \item{New operations companies: }{Average spread of new credit operations - Companies - Total - p.p.}
#'  \item{New operations households: }{Average spread of new credit operations - Households - Total - p.p.}
#'  \item{Non-revolving total: }{Average spread of new non-revolving credit operations - Total - p.p.}
#'  \item{Non-revolving companies: }{Average spread of new non-revolving credit operations - Companies - Total - p.p.}
#'  \item{Non-revolving households: }{Average spread of new non-revolving credit operations - Households - Total - p.p.}
#'}
#' @details For Indicator = 'House debt', the parameter Index can be Total and Without mortgage.
#'\itemize{
#'  \item{Total: }{Household debt - \%.}
#'  \item{Without mortgage: }{Household debt without mortgage loans - \%.}
#'}
#'
#' @details For Indicator = 'Interest rate', the parameter Index can be Total, New operations companies, New operations households, Non-revolving total, Non-revolving companies and Non-revolving households.
#'
#'\itemize{
#'  \item{Total: }{Average interest rate of new credit operations - Total - \% (per year) p.y.}
#'  \item{New operations companies: }{Average interest rate of new credit operations - Companies - Total - \% p.y.}
#'  \item{New operations households: }{Average interest rate of new credit operations - Households - Total - \% p.y.}
#'  \item{Non-revolving total: }{Average interest rate of new non-revolving credit operations - Total - \% p.y.}
#'  \item{Non-revolving companies: }{Average interest rate of new non-revolving credit operations - Companies - Total - \% p.y.}
#'  \item{Non-revolving households: }{Average interest rate of new non-revolving credit operations - Households - Total - \% p.y.}
#'}
#'
#' @return A list of objects containing the dataset with the time series associated to the Indicator and Index, a vector with the dates of the times series, the Indicator and the Index.
#' @export
#' @importFrom reshape2 "melt"
#' @importFrom dplyr "select" "filter" "%>%"
#' @importFrom utils "head" "read.csv"
#' @seealso \code{\link{FitModels}} and \code{\link{plot_FitModels}}.
#' @examples
#' data0 <- ReadData(Indicator = 'Interest rate', Index = 'Total')
#' data1 <- ReadData(Indicator = 'Spread', Index = 'New operations companies')
#' data2 <- ReadData(Indicator = 'Household debt', Index = 'Without mortgage')
#' data3 <- ReadData(Indicator = 'Allowance', Index = 'Public banks')

# library(dplyr)
# library(reshape2)

ReadData <- function(Indicator, Index) {

  Aux_files <- c('Allowance_198806_201903.csv',
                 '90_past_due_loans_201103_201903.csv',
                 '90_past_due_loans_by_origin_capital_20003_201903.csv',
                 'Average_spread_yearly_201103_201903.csv',
                 'House_debt_200501_201903.csv',
                 'Interest_rates_yearly_201103_201903.csv')

  Aux_info <- c('Allowance',
                'Due loans',
                'Due loans by capital',
                'Spread',
                'Household debt',
                'Interest rate')
  Index_aux = ''

  if (Indicator == 'Allowance' && Index=='Total')                   {Index_aux = 'X13645'} else
  if (Indicator == 'Allowance' && Index=='Public banks')            {Index_aux = 'X13666'} else
  if (Indicator == 'Allowance' && Index=='National private control'){Index_aux = 'X13672'} else
  if (Indicator == 'Allowance' && Index=='Foreign control')         {Index_aux = 'X13678'} else
  if (Indicator == 'Allowance' && Index=='Private control')         {Index_aux = 'X13684'}

  if (Indicator == 'Spread' && Index=='Total')                      {Index_aux = 'X20783'} else
  if (Indicator == 'Spread' && Index=='New operations companies')   {Index_aux = 'X20784'} else
  if (Indicator == 'Spread' && Index=='New operations households')  {Index_aux = 'X20785'} else
  if (Indicator == 'Spread' && Index=='Non-revolving total')        {Index_aux = 'X27631'} else
  if (Indicator == 'Spread' && Index=='Non-revolving companies')    {Index_aux = 'X27632'} else
  if (Indicator == 'Spread' && Index=='Non-revolving households')   {Index_aux = 'X27633'}

  if (Indicator == 'Household debt' && Index=='Total')                  {Index_aux = 'X20400'} else
  if (Indicator == 'Household debt' && Index=='Without mortgage')       {Index_aux = 'X19882'}

  if (Indicator == 'Interest rate' && Index=='Total')                     {Index_aux = 'X20714'} else
  if (Indicator == 'Interest rate' && Index=='New operations companies')  {Index_aux = 'X20715'} else
  if (Indicator == 'Interest rate' && Index=='New operations households') {Index_aux = 'X20716'} else
  if (Indicator == 'Interest rate' && Index=='Non-revolving total')       {Index_aux = 'X27623'} else
  if (Indicator == 'Interest rate' && Index=='Non-revolving companies')   {Index_aux = 'X27624'} else
  if (Indicator == 'Interest rate' && Index=='Non-revolving households')  {Index_aux = 'X27625'}

  if (Index_aux == '') {stop("ERROR: At least one of the parameters (Indicator or Index) was misspelled. Check it out and try again.")}

  aux <- as.data.frame(cbind(Aux_files, Aux_info))

  # Selecting the databaset associated to the indicator
  NamesFiles <- aux %>% filter(Aux_info==Indicator) %>% select(Aux_files)

  # Reading the data and deleting the last row
  db <- read.csv(text=paste0(head(readLines(system.file("extdata", paste('',NamesFiles[1,1],'', sep =''), package = "CredInd")), -1), collapse="\n"), sep = ';', header=T, na.strings = c('-'))

  # Changing the names of the columns
  names(db) <- c(substr(names(db),1,6))

  # Transposing the data
  db <- db %>%
        melt(id.vars = c('Date'), na.rm = T) %>%
        filter(variable == Index_aux)

  db$variable <- paste(Indicator,'-', Index)

  final_list <- list(data = db,
                     Date = db$Date,
                     Indicator = Indicator,
                     Index = Index)

  class(final_list) = 'CredInd'

  return(final_list)
}

# ReadData(Indicator = 'bababa rate', Index = 'Total')
# data1 <- ReadData(Indicator = 'Interest rate', Index = 'Total')
# ReadData(Indicator = 'Interest rate', Index = 'New operations companies')
# ReadData(Indicator = 'Interest rate', Index = 'New operations households')
# ReadData(Indicator = 'Interest rate', Index = 'Non-revolving total')
# ReadData(Indicator = 'Interest rate', Index = 'Non-revolving companies')
# ReadData(Indicator = 'Interest rate', Index = 'Non-revolving households')
#
# ReadData(Indicator = 'Spread', Index = 'Total')
# ReadData(Indicator = 'Spread', Index = 'New operations companies')
# ReadData(Indicator = 'Spread', Index = 'New operations households')
# ReadData(Indicator = 'Spread', Index = 'Non-revolving total')
# ReadData(Indicator = 'Spread', Index = 'Non-revolving companies')
# ReadData(Indicator = 'Spread', Index = 'Non-revolving households')
#
# ReadData(Indicator = 'Household debt', Index = 'Total')
# ReadData(Indicator = 'Household debt', Index = 'Without mortgage')
#
# ReadData(Indicator = 'Allowance', Index = 'Total')
# ReadData(Indicator = 'Allowance', Index = 'Public banks')
# ReadData(Indicator = 'Allowance', Index = 'National private control')
# ReadData(Indicator = 'Allowance', Index = 'Foreign control')
# ReadData(Indicator = 'Allowance', Index = 'Private control')
