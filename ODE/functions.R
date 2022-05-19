set.seed(1)
library(ZVCV)

y <- read.csv("y.csv", header = FALSE) ## observed data
colnames(y) <- NULL
rownames(y) <- NULL
y <- as.numeric(y)
n <- length(y) ## n = 11

F <- function(x, theta){
  return(c(x[2], (theta*(1 - x[1]^2)*x[2]) - x[1]))
}

## Runge Kutta method for solving ODE
RK <- function(theta){
  n_int <- 100
  x <- matrix(0, nrow = 2, ncol = 1+n_int)
  x[,1] <- c(0,2)
  h <- 10/n_int
  for (i in 1:n_int){
    k1 <- F(x[,i], theta)
    k2 <- F(x[,i] + (h*k1/2), theta)
    k3 <- F(x[,i] + (h*k2/2), theta)
    k4 <- F(x[,i] + h*k3, theta)
    x[,(i+1)] <- x[,i] + h*(k1 + 2*k2 + 2*k3 + k4)/6
  }
  iter <- 1 + seq(0, 10, 1)*n_int/10
  x_a <- x[1, iter]
  return(x_a)
}



## Calculates d x_a/d theta
RK_deriv <- function(theta, k){
  return((RK(theta + k) - RK(theta))/k)
}


## Calculates log likelihood: p(y|theta)
loglikelihood<- function(y, theta, sigma, x_a){
  prod <- 0
  n <- length(y)
  for (j in 1:n){
    prod <- prod + log(dnorm(as.numeric(y[j]), mean = x_a[j], sd = sigma))
  }
  prod <- prod - (n/2)*log(2*pi*sigma^2)
  return(prod)
}

## Calculates power density: p(theta|y,t) = p(y|theta)^t p(theta)
power_density <- function(theta, y, t, sigma, x_a){
  return((exp(loglikelihood(y, theta, sigma, x_a))^t)*dlnorm(theta, meanlog = 0, sdlog = 0.5))
}

## Calculates gradient of power density 
power_grad <- function(theta, k, t, y,sigma, x_a){
  n <- length(y)
  x_a_prime <- RK_deriv(theta, k)
  u <- crossprod(y - x_a, x_a_prime)
  u <- (t/sigma^2)*u
  u <- u - ((1 + 4*log(theta))/theta)
  return(u)
}

## Metropolis-Hastings algorithm for sampling fromm  power density
## with proposal distribution equal to prior distribution
MH <- function(t, N, sigma, y){
  thetas <- rep(0, N)
  accept <- 1
  itr <- 0
  thetas[1] <- rnorm(1, 1, 0.05)
  while(accept < N){
    itr <- itr+1
    prop <- rlnorm(1, meanlog = log(thetas[accept]), sdlog = 0.5)
    x_a <- RK(prop)
    ratio <- power_density(prop, y, t, sigma, RK(prop))/power_density(thetas[accept], y, t, sigma, RK(thetas[accept]))
    if (!is.na(ratio))
      ratio <- min(1, ratio)
    else
      ratio <- 0
    #print(ratio)
    u <- runif(1)
    if(ratio >= u){
      accept = accept + 1
      thetas[accept] = prop
    }
  }
  print(accept/itr)
  return(thetas)
}

