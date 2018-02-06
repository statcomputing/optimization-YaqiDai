
## function to calculate the MLE of Cauchy distribution by Newton-Raphson method
NewtonMleCauchy <- function(X, a, tolar){
  # X is vector of observed sample
  # a is the start value
  n <- length(X)
  theta.iterate <- a
  #calculate the first derivate of log likelyhood of Cauchy distribution
  firstd <- (-2)*sum((theta.iterate - X)/(1 + (theta.iterate - X)^2))
  while(abs(firstd) > tolar ){
    secondd <- (-2)*sum(((1 - (theta.iterate - X)^2)/(1 + (theta.iterate - X)^2))^2)
    theta.new <- theta.iterate - firstd/secondd
    theta.iterate <-theta.new
    firstd <- (-2)*sum((theta.iterate - X)/(1 + (theta.iterate - X)^2))
  }
  return(theta.iterate)
}

X <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
       3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
Y <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
theta <- matrix(NA, nrow = length(Y), ncol = 1)
for(i in 1:length(Y)){
  theta[i] <- NewtonMleCauchy(X, Y[i], 1E-10)
  }

