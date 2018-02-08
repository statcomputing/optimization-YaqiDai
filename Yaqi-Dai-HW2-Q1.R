################################################
###(b)
## log-likelihood function

Loglikelihood <- function(X, theta){
  n <- length(X)
  -n*log(pi) - sum(log(1 + (theta - X)^2))
}

## function to calculate the MLE of Cauchy distribution by Newton-Raphson method
NewtonMleCauchy <- function(X, a, tolar){
  # X is vector of observed sample
  # a is the start value
  timestart <- Sys.time()
  n <- length(X)
  theta.iterate <- a
  i <- 0 # i is the number of iteration
  #calculate the first derivate of log likelyhood of Cauchy distribution
  firstd <- (-2)*sum((theta.iterate - X)/(1 + (theta.iterate - X)^2))
  while(abs(firstd) > tolar ){
    secondd <- (-2)*sum(((1 - (theta.iterate - X)^2)/(1 + (theta.iterate - X)^2))^2)
    theta.new <- theta.iterate - firstd/secondd
    theta.iterate <-theta.new
    firstd <- (-2)*sum((theta.iterate - X)/(1 + (theta.iterate - X)^2))
    i <- i + 1
  }
  timeend <- Sys.time()
  running.time <- timeend - timestart
  result1 <- c(theta.iterate, running.time, i)
  return(result1)
}

X <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
       3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
Y <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)

## plot log-lilihood
theta <- seq(-20, 20, by = 0.01)
loglikelihood <- sapply(theta, FUN=function(theta) Loglikelihood(X, theta))
plot(theta, loglikelihood, type="l")

## find MLE for theta
mle.metrix <- sapply(Y, FUN = function(Y) NewtonMleCauchy(X, Y, 1E-10))

## start at sample mean
sample.mean <- mean(X)
theta.sm <- NewtonMleCauchy(X, sample.mean, 1E-10)

#############################################################################
### (c)

## function for Fixedpoint
Fixedpoint <- function(alfa, a, totlar){
  # a is the initial value of theta
  timestart <- Sys.time()
  theta.iterate <- a
  firstd <- (-2)*sum((theta.iterate - X)/(1 + (theta.iterate - X)^2))
  theta.new <- theta.iterate + alfa*firstd
  i <- 0
  while (abs(firstd) > totlar && i <= 2000) {
    theta.new <- theta.iterate + alfa*firstd
    theta.iterate <- theta.new
    firstd <- (-2)*sum((theta.iterate - X)/(1 + (theta.iterate - X)^2))
    i <- i + 1
  }
  timeend <- Sys.time()
  runningtime <- timestart - timeend
  result2 <- c(theta.iterate, runningtime, i)
  return(result2)
}

alfa.set <- c(1, 0.64, 0.25)

theta1 <- sapply(Y, FUN = function(Y)Fixedpoint(1, Y, 1E-10))
theta2 <- sapply(Y, FUN = function(Y)Fixedpoint(0.64, Y, 1E-10))
theta3 <- sapply(Y, FUN = function(Y)Fixedpoint(0.25, Y, 1E-10))

##########################################################################
###(d)

## function to calculate the MLE of Cauchy distribution by Fisher scoring
FisherScoring <- function(X, a, tolar){
  # X is vector of observed sample
  # a is the start value
  n <- length(X)
  fisherinfor <- n/2
  theta.iterate <- a
  #calculate the first derivate of log likelyhood of Cauchy distribution
  firstd <- (-2)*sum((theta.iterate - X)/(1 + (theta.iterate - X)^2))
  i <- 0
  while(abs(firstd) > tolar && i <= 20){
    theta.new <- theta.iterate + firstd/fisherinfor
    theta.iterate <-theta.new
    firstd <- (-2)*sum((theta.iterate - X)/(1 + (theta.iterate - X)^2))
    i <- i + 1
  }
  return(theta.iterate)
}
## use theta calculated from Fisher-scoring method as the initial value 
## in the Newton-Raphson method

NewtonFisher <- function(X, a, tolar){
  theta.fisher <- FisherScoring(X, a, tolar)
  theta.newton <- NewtonMleCauchy(X, theta.fisher, tolar)
  return(theta.newton)
}

## theta estimator at different initial value


fisher.estimate <- matrix(sapply(Y, FUN = function(Y)FisherScoring(X, Y,tolar = 1E-10)))
newton.estimate <- sapply(Y, FUN = function(Y)NewtonFisher(X, Y, tolar = 1E-10))
