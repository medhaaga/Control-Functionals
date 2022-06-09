dist <- function(x1, x2, theta1, theta2){
  ret <- (1 + theta1*sum(x1^2))*(1 + theta1*sum(x2^2))
  ret <- exp(-(1/2*(theta2^2))*sum((x1-x2)^2))/ret
  return(ret)
}

ones <- function(n){
  return(matrix(1, nrow = n, ncol = 1))
}

vecA <- function(data, m, theta1, theta2){
  n <- nrow(data)
  K0 <- matrix(0, m, m)
  for (i in 1:m){
    for (j in 1:i){
      K0[i,j] <- dist(data[i,], data[j,], theta1, theta2)
      K0[j,i] <- K0[i,j]
    }
  }
  K0 <- solve(K0 + lam*diag(m))
  
  K1 <- matrix(0, n-m, m)
  for (i in 1:(n-m)){
    for (j in 1:m){
      K1[i,j] <- dist(data[(m+i),], data[j,], theta1, theta2)
    }
  }
  
  K_mult <- K1 %*% K0
  
  term1 <- (as.numeric(t(ones(n-m))%*%K_mult%*%ones(m))/(1+as.numeric(t(ones(m))%*%K0%*%ones(m))))*(K0%*%ones(m))
  term2 <- t(K_mult)%*%ones(n-m)
  A <- term1 - term2
  
  return(A)
  
}

log.density <- function(x, p, mu1, mu2, sd1, sd2){
  return(log(p*dnorm(x, mean = mu1, sd = sd1) + (1-p)*dnorm(x, mean = mu2, sd= sd2)))
}


mh.mcmc <- function(start, p, mu1, mu2, sd1, sd2, N, h){
  X <- rep(0, N)
  X[1] <- start
  for (i in 2:N){
    prop <- rnorm(1, mean = X[i-1], sd = 2)
    ratio <- log.density(prop, p, mu1, mu2, sd1, sd2) - log.density(X[i-1], p, mu1, mu2, sd1, sd2)
    if(runif(1) < exp(ratio))
      X[i] <- prop else
        X[i] <- X[i-1]
  }
  return (X)
}

iid <- function(start, p, mu1, mu2, sd1, sd2, N, h){
  X <- rep(0, N)
  
  for (i in 1:N){
    bern <- rbinom(1, size=1, p <- p)
    if(bern == 1) X[i] <- rnorm(1, mean = mu1, sd = sd1)
    else X[i] <- rnorm(1, mean = mu2, sd = sd2)
  }
  return (X)
}


u <- function(x, p, mu1, mu2, sd1, sd2){
  dens <- exp(log.density(x, p, mu1, mu2, sd1, sd2))
  grad <- -((p*(x - mu1)*(1/sd1^2)*dnorm(x, mu1, sd1)) + ((1-p)*(x - mu2)*(1/sd2^2)*dnorm(x, mu2, sd2)))
  return(grad/dens)
}

rep <- 100

mu_AM <- rep(0, rep)
mu_CF <- rep(0, rep)
mu_sCF <- rep(0, rep)
mu_ICF <- rep(0, rep)
mu_ZVCV <- rep(0, rep)
p <- 0.5
mu1 <- -2
mu2 <- 2
sd1 <- 1
sd2 <- 1
N <- 100
m <- N/2
theta1 <- 0.1
theta2 <- 1

for (i in 1:rep){
  if(i%%10 == 0) print(i)
  
  is.samp <- rnorm(N, sd = 2)
  weights <- (p*dnorm(is.samp, mean = mu1, sd=  sd1) + (1-p)*dnorm(is.samp, mean = mu2, sd = sd2))/dnorm(is.samp, sd=2)
  grad <- as.matrix(apply(as.matrix(is.samp, N, 1), 1, u, p, mu1, mu2, sd1, sd2), nrow = N, ncol = 1)
  integrands <- as.matrix(exp(is.samp), nrow = N, ncol = 1)
  
  mu_AM[i] <- sum((weights/sum(weights))*integrands)
  
  A <- vecA(as.matrix(is.samp, nrow = N, ncol = 1), m, theta1, theta2)
  W0 <- diag(weights[1:m])/sum(weights[1:m])
  W1 <- diag(weights[(m+1):N])/sum(weights[(m+1):N])
  mu_ICF[i] <- t(ones(N-m))%*%W1%*%integrands[(m+1):N,] + (1/(N-m))*t(A)%*%W0%*%integrands[1:m,]
}
boxplot(data.frame("AM" = mu_AM, "ICF" = mu_ICF))



