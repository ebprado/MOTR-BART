library(CredInd)

context('Testing the possible inputs and outputs for all the functions of the package test')

test_that('Testing the function RealData', {

  expect_error(ReadData(Indicator = 'bla blab rate', Index = 'Total')) # "ERROR: At least one of the parameters (Indicator or Index) was misspelled. Check it out and try again."
  expect_error(ReadData(Indicator = 'Interest rate', Index = 'bla b')) # "ERROR: At least one of the parameters (Indicator or Index) was misspelled. Check it out and try again."
  expect_output(str(ReadData(Indicator = 'Interest rate', Index = 'Total')), "List of 4")

})

test_that('Testing the function FitModels', {

  data1 <- ReadData(Indicator = 'Interest rate', Index = 'Total')
  expect_output(str(FitModels(data=data1, horizon=20)), "List of 6")
  expect_error(FitModels(data=data1, horizon='20'), "ERROR: Either the class of the data is not 'CredInd' or the argument horizon is not numeric!")
  expect_error(FitModels(data=rnorm(100), horizon=20), "ERROR: The class of the data must be 'CredInd'. Generate the data by using the function ReadData")
  expect_error(FitModels(data=rnorm(100), horizon='20'), "ERROR: Either the class of the data is not 'CredInd' or the argument horizon is not numeric!")

})

test_that('Testing the function plot_FitModels', {

  data1 <- ReadData(Indicator = 'Interest rate', Index = 'Total')
  mod1 <- FitModels(data=data1, horizon=20)
  expect_output(str(plot_FitModels(mod1)), NULL)
  expect_error(plot_FitModels('mod1')) # "ERROR: The class of the data is not 'CredInd'. Generate the data by using the function FitModels."

})
