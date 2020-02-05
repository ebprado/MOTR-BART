# This script provides a brief introduction to:
# 1. Principal Components Analysis
# 2. Discriminant Analysis
# 3. Hierarchical clustering
# 4. Model-based clustering 

##########################################################################
##########################################################################
# PRINCIPAL COMPONENTS ANALYSIS (PCA)
##########################################################################
##########################################################################

require(graphics)
data(mtcars)
head(mtcars)
library(klaR)
# [, 1]	mpg	Miles/(US) gallon
# [, 2]	cyl	Number of cylinders
# [, 3]	disp	Displacement (cu.in.)
# [, 4]	hp	Gross horsepower
# [, 5]	drat	Rear axle ratio
# [, 6]	wt	Weight (1000 lbs)
# [, 7]	qsec	1/4 mile time
# [, 8]	vs	Engine (0 = V-shaped, 1 = straight)
# [, 9]	am	Transmission (0 = automatic, 1 = manual)
# [,10]	gear	Number of forward gears
# [,11]	carb	Number of carburetors

boxplot(mpg~cyl, data=mtcars, ylab='Miles per gallon (mpg)', xlab='Number of cylinders', main= 'Cylinders x Mpg')
boxplot(mpg~carb, data=mtcars, ylab='Miles per gallon (mpg)', xlab='Number of carburetors' , main= 'Carburetors x Mpg')
boxplot(hp~carb, data=mtcars, ylab='Horse power', xlab='Number of carburetors' , main= 'HP x Carburetors')
boxplot(mpg~gear, data=mtcars, ylab='Miles per gallon (mpg)', xlab='Number of gears' , main= 'Gears x Mpg')
plot(mpg~hp, data=mtcars, ylab='Miles per gallon (mpg)', xlab='Horse power (hp)' , main= 'HP x Mpg')

## Running Principal Components Analysis (PCA)

pca <- prcomp(~ mpg + cyl + disp + hp + drat + wt + qsec + vs + am + gear + carb, data = mtcars, scale=FALSE)
pca
summary(pca)

## Scale = TRUE standardises the values - important if measured in different units  

pca <- prcomp(~ mpg + cyl + disp + hp + drat + wt + qsec + vs + am + gear + carb, data = mtcars, scale=TRUE)
pca
summary(pca)
biplot(pca)
pca$scale

## Creating the PCs
scale_mtcars <- as.matrix(scale(mtcars))
pcs <- scale_mtcars%*%as.matrix(pca$rotation)
cor(pcs)

##########################################################################
##########################################################################
# Linear Discriminant Analysis (LDA) #####################################
##########################################################################
##########################################################################

# An idea about LDA
##########################################################################
n <- 1000
red <- c(rnorm(n*0.6, -2, 1/4), rnorm(n*0.4, 0, 1))
blue <-c(rnorm(n*0.6, 0, 1/4), rnorm(n*0.4, 2, 1))

db1 <- data.frame(label = rep('red', n), y = rnorm(n, 4, 1), x = red)
db2 <- data.frame(label = rep('blue', n), y = rnorm(n, 0, 1), x = blue)
db <- rbind(db1, db2)

# aa <- lda(label ~ ., db, CV = TRUE)
aa <- lda(label ~ y + x, db, CV = FALSE)
aa_values <- predict(aa)
# aa_values$x
# plot(aa_values$posterior[,2], aa_values$x)
aa$scaling
aa$posterior
aux <- rbind(db1, db2)
ggplot(aux, aes(x=x, y=y, color=label)) +
  geom_point() 
# geom_abline(intercept = 2, slope=0.3316098)
partimat(label~.,data=db,method="lda") 

##########################################################################

# Use the Iris data for this example
data(iris)

# Have a look at the data set
str(iris)
# So there are four measurements on each of a set of different plants

# The plants are:
table(iris$Species)

# The idea of LDA is that it builds a classification model based on a straight line of best fit that separates out the different groups
# Here we can use LDA to help us predict the Species variable from the petal and sepal measurements

# Suppose we only knew 100 of the measurements
iris_train_rows = sample(1:nrow(iris), 100)

# So our data set was really
iris_train = iris[iris_train_rows,]

# We build an LDA classifier with
lda1 = lda(Species ~ ., data = iris_train)

# Look at the output
# print(lda1)
partimat(Species~.,data=iris,method="lda") 
# Like the PC above this creates straight lines of discrimination (here they're called 'coefficients of linear discrimination')
# These allow us to classify the three groups depending on which side of the line they're on

# You can visualise the performance with:
plot(lda1) 
# You should it separate oing out some of the values
# You can also see that setosa is well separated (and thus easily identifiable) whilst the others are harder to separate

# For the left out data we can create predictions and evaluate how well it did

# First create the test data
iris_test_rows = (1:nrow(iris))[-iris_train_rows]
iris_test = iris[iris_test_rows,]

# Now create the predictions
predicted = predict(lda1, iris_test)

# Now tabulate the true values against the predicted
table(iris_test$Species, predicted$class)
# It got all of them right apart from 2!


##########################################################################
##########################################################################
# Hierarchical clustering (hclust) ---------------------------------------
##########################################################################
##########################################################################

# Go back to the arrests data
# The idea of clustering is that you group together rows of the data set
# For this example that means finding states which are similar to each other

# This method is called hierarchical clustering because it starts with each state in a separate cluster, 
# then slowly joins them together one by one until they are all in the same cluster

# Reminder of the data set
str(USArrests)

# The first step is to calculate the distance between each state. We can do this using the dist function
my_dist = dist(USArrests)

# If you print this out you get a big matrix which shows the 'distance' between each state
# (Note that distance here means distance between murder/rape/urbanPop/assault values not location distance)
round(my_dist,2)

# You now run the hclust on this
hc1 = hclust(my_dist)

# Plot it
plot(hc1)
# This is a dendrogram
# It shows you that there are two broad clusters. You can then try to interpret why each one is put together. 
# People often choose the number of clusters based on the 'persistence' of the dendrogram. Here for example
# two clusters seems to take up most of the distance on the dendrogram

##########################################################################
##########################################################################
# Model-based clustering (Mclust) ----------------------------------------
##########################################################################
##########################################################################

# Mclust is a much fancier version of clustering that is model-based, i.e. has a probability model 
# underneath it so that we can give a probability to an observation being in a certain cluster
# It is not hierarchical so you need to choose the number of clusters before you run it
# Often you run it for a range of clusters and Mclust returns the BIC to help you choose
# the 'optimum' number of clusters

# You can do this in one go for e.g. the USArrests data
mc1 = Mclust(USArrests)

# Summarise it
summary(mc1) 
# So it prefers using 3 clusters for this dat set
# The first cluster has 20 states in it, the second 10, and the third 20

# Plot it
plot(mc1,  what = "BIC")
# This shows the number of clusters and the (negative) BIC. The 'best' model has the highest BIC
# The different lines show different types of covariance matrices used for the model
# (You can probably ignore them!)

# This shows the different cluster for each of the variables
plot(mc1, what = "classification")

# You could also try doing this for the Iris data

# Final job! --------------------------------------------------------------

# Now get these working on your own data sets!
