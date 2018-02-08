#####################################
###(a)
beetles <- data.frame(
  days <- c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154),
  nbeetles <- c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024)
)
K <- 1500
r1 <- log(beetles$nbeetles*(K - 2)/(K - beetles$nbeetles)*2)
r2 <- r1 / beetles$days....c.0..8..28..41..63..69..97..117..135..154.
r.mean <- mean(r2[2: 10])
N0 <- beetles[1,2]

nls(nbeetles ~ (K*N0)/(N0 + (K - N0)*exp(-r*days)),
    start = list(K = 1500, r = r.mean), data = beetles, trace = TRUE)


##################################################
####(b)

K1 <- seq(500, 1500, length.out = 200)
r1 <- seq(0.1, 0.2, length.out = 200)

n.beetles <- as.numeric(beetles$nbeetles....c.2..47..192..256..768..896..1120..896..1184..1024.)
n.days <- as.numeric(beetles$days....c.0..8..28..41..63..69..97..117..135..154.)

SSE <- function(K1, r1){
  error.sq <- sum((n.beetles - (N0*K1)/(N0 + (K1 - N0)*exp(-r1*n.days)))^2)
  return(error.sq)
}

matrix.a <- matrix(nrow = 200, ncol = 200)

for(i in 1:200){
  for(j in 1:200){
    matrix.a[i, j] <- SSE(K1[i], r1[j])
  }
}

filled.contour(K1, r1, matrix.a, plot.title = title(main = "Contour plot of the sum squared errors", 
                                                    xlab = "K", ylab = "r"),
               key.title = title(main = "SSE"))
###########################################################
###(C)
nloglikelihood.c <- function(theta, N, days){
  K <- theta[1]
  r <- theta[2]
  sigma <- theta[3]
  t <- days
  mu <- log((K*N0)/(2 + (K - N0)*exp(-r*t)))
  theta.iterate <- (- sum(dnorm(log(N), mu, sigma, log = TRUE)))
}


sq <- sqrt(var(log(beetles$nbeetles)))

theta0 <- c(K, r.mean, sq)
mle <- nlm(nloglikelihood.c, theta0, N = beetles$nbeetles,days = beetles$days,hessian = TRUE)
estimate <- mle$estimate
estimate
matrix.hessian <- mle$hessian
matrix.hessian
variance <- solve(matrix.hessian)
variance
diag(variance)
