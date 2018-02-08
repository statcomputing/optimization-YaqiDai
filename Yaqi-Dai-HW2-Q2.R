X <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
       2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)
#######################################################################
## Q2 (a)

## log-likelihood function
Loglikelihood <- function(X, theta){
  n <- length(X)
  sum(log(1 - cos(X - theta))) - n*log(2*pi)
}

X <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
       2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)
theta <- seq(-pi, pi, by = 0.01)

plot(theta, sapply(theta, FUN=function(theta) Loglikelihood(X, theta)), type="l")

################################################################################
## Q2 (b)
xmean <- mean(X)
theta.moment <- asin(xmean - pi)

###################################################################################
## Q2(c)

NewtonMethod <- function(X, a, tolar){
  # X is vector of observed sample
  # a is the start value
  n <- length(X)
  theta.iterate <- a
  #calculate the first derivate of log likelyhood of Cauchy distribution
  firstd <- sum(sin(X - theta.iterate)/(cos(X - theta.iterate) - 1))
  while(abs(firstd) > tolar ){
    secondd <- sum(1/(cos(X - theta.iterate) - 1))
    theta.new <- theta.iterate - firstd/secondd
    theta.iterate <-theta.new
    firstd <- sum(sin(X - theta.iterate)/(cos(X - theta.iterate) - 1))
  }
  return(theta.iterate)
}

theta.c <- NewtonMethod(X, theta.moment, 1E-10)
## Q2(d)

theta.d1 <- NewtonMethod(X, -2.7, 1E-10)
theta.d2 <- NewtonMethod(X, 2.7, 1e-10)

## Q2(e)

initialvalue <- seq(-pi, pi, length.out = 200)
theta.set <- round(sapply(initialvalue, 
                          FUN = function(initialvalue)
                            NewtonMethod(X, initialvalue, 1E-10)), 4)
theta.list <- split(initialvalue, theta.set)
list(theta.list)
