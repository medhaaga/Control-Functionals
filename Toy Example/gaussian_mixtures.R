library(ZVCV)
set.seed(1)
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
mu_ZVCV <- rep(0, rep)
N <- 200

###########################################################
###### Case 1: f(x) = exp(x) p(x) = N(0,1) ################
###########################################################

p <- 0
mu1 <- 0
mu2 <- 0
sd1 <- 1
sd2 <- 1

for (i in 1:rep){
  if(i%%10 == 0) print(i)
  samples <- as.matrix(mh.mcmc(start = -3, p, mu1, mu2, sd1, sd2, N+500, 1)[501:(N+500)], nrow = N, ncol = 1)
  grad <- as.matrix(apply(samples, 1, u, p, mu1, mu2, sd1, sd2), nrow = n_samp, ncol = 1)
  integrands <- as.matrix(exp(samples), nrow = n_samp, ncol = 1)
  mu_AM[i] <- mean(integrands)
  mu_ZVCV[i] <- zvcv(integrand = integrands, samples = samples, derivatives = grad)$expectation
  
  samples <- unique(samples)
  n_samp <- nrow(samples)
  grad <- as.matrix(apply(samples, 1, u, p, mu1, mu2, sd1, sd2), nrow = n_samp, ncol = 1)
  integrands <- as.matrix(exp(samples), nrow = n_samp, ncol = 1)
  
  mu_sCF[i] <- CF(integrands = integrands, samples = samples, derivatives = grad, steinOrder = 1, kernel_function = "prodsim", sigma = c(.1,1))$expectation
  ind <- sample.int(n_samp, n_samp/2, replace = FALSE)
  cf <- CF(integrands = integrands, samples = samples, derivatives = grad, steinOrder = 1, kernel_function = "prodsim", sigma = c(0.1,1), est_inds = ind)
  mu_CF[i] <- cf$expectation
  
  mu_ZVCV[i] <- zvcv(integrand = integrands, samples = samples, derivatives = grad)$expectation
}

pdf(file = "gaussian_mixture_mix0_exp.pdf", height = 5, width = 5)
boxplot(data.frame("AM" = mu_AM, "ZV-CV" = mu_ZVCV, "simplified CF" = mu_sCF), ylab = expression(mu(f)))
dev.off()

truth <- exp(0.5)
mse <- round(colMeans((cbind(mu_AM, mu_ZVCV, mu_sCF) - truth)^2),4)
print(paste(round(mean(mu_AM),4), mse[1], round(mean(mu_ZVCV), 4), mse[2], round(mean(mu_CF), 4), mse[3]))

###########################################################
### Case 2: f(x) = x^2 p(x) = 0.5 N(-1,1) + 0.5 N(1,1) ####
###########################################################

p <- 0.5
mu1 <- -1
mu2 <- 1
sd1 <- 1
sd2 <- 1

for (i in 1:rep){
  if(i%%10 == 0) print(i)
  samples <- as.matrix(mh.mcmc(start = -3, p, mu1, mu2, sd1, sd2, N+500, 1)[501:(N+500)], nrow = N, ncol = 1)
  grad <- as.matrix(apply(samples, 1, u, p, mu1, mu2, sd1, sd2), nrow = n_samp, ncol = 1)
  integrands <- as.matrix((samples^2), nrow = n_samp, ncol = 1)
  mu_AM[i] <- mean(integrands)
  mu_ZVCV[i] <- zvcv(integrand = integrands, samples = samples, derivatives = grad)$expectation
  
  samples <- unique(samples)
  n_samp <- nrow(samples)
  grad <- as.matrix(apply(samples, 1, u, p, mu1, mu2, sd1, sd2), nrow = n_samp, ncol = 1)
  integrands <- as.matrix((samples^2), nrow = n_samp, ncol = 1)
  
  mu_sCF[i] <- CF(integrands = integrands, samples = samples, derivatives = grad, steinOrder = 1, kernel_function = "prodsim", sigma = c(.1,1))$expectation
  ind <- sample.int(n_samp, n_samp/2, replace = FALSE)
  mu_CF[i] <- CF(integrands = integrands, samples = samples, derivatives = grad, steinOrder = 1, kernel_function = "prodsim", sigma = c(0.1,1), est_inds = ind)$expectation
  mu_ZVCV[i] <- zvcv(integrand = integrands, samples = samples, derivatives = grad)$expectation
}

pdf(file = "gaussian_mixture_mix1_poly.pdf", height = 5, width = 5)
boxplot(data.frame("AM" = mu_AM, "ZV-CV" = mu_ZVCV, "simplified CF" = mu_sCF), ylab = expression(mu(f)))
dev.off()

truth <- 2
mse <- round(colMeans((cbind(mu_AM, mu_ZVCV, mu_sCF) - truth)^2),4)
print(paste(round(mean(mu_AM),4), mse[1], round(mean(mu_ZVCV), 4), mse[2], round(mean(mu_CF), 4), mse[3]))

###########################################################
### Case 3: f(x) = x^2 p(x) = 0.5 N(-2,1) + 0.5 N(2,1) ####
###########################################################

p <- 0.5
mu1 <- -2
mu2 <- 2
sd1 <- 1
sd2 <- 1

for (i in 1:rep){
  if(i%%10 == 0) print(i)
  samples <- as.matrix(mh.mcmc(start = -3, p, mu1, mu2, sd1, sd2, N+500, 1)[501:(N+500)], nrow = N, ncol = 1)
  grad <- as.matrix(apply(samples, 1, u, p, mu1, mu2, sd1, sd2), nrow = n_samp, ncol = 1)
  integrands <- as.matrix((samples^2), nrow = n_samp, ncol = 1)
  mu_AM[i] <- mean(integrands)
  mu_ZVCV[i] <- zvcv(integrand = integrands, samples = samples, derivatives = grad)$expectation
  
  samples <- unique(samples)
  n_samp <- nrow(samples)
  grad <- as.matrix(apply(samples, 1, u, p, mu1, mu2, sd1, sd2), nrow = n_samp, ncol = 1)
  integrands <- as.matrix((samples^2), nrow = n_samp, ncol = 1)
  
  mu_sCF[i] <- CF(integrands = integrands, samples = samples, derivatives = grad, steinOrder = 1, kernel_function = "prodsim", sigma = c(.1,1))$expectation
  ind <- sample.int(n_samp, n_samp/2, replace = FALSE)
  mu_CF[i] <- CF(integrands = integrands, samples = samples, derivatives = grad, steinOrder = 1, kernel_function = "prodsim", sigma = c(0.1,1), est_inds = ind)$expectation
  mu_ZVCV[i] <- zvcv(integrand = integrands, samples = samples, derivatives = grad)$expectation
}

print()

pdf(file = "gaussian_mixture_mix2_poly.pdf", height = 5, width = 5)
boxplot(data.frame("AM" = mu_AM, "ZV-CV" = mu_ZVCV, "simplified CF" = mu_sCF), ylab = expression(mu(f)))
dev.off()

truth <- 5
mse <- round(colMeans((cbind(mu_AM, mu_ZVCV, mu_sCF) - truth)^2),4)
print(paste(round(mean(mu_AM),4), mse[1], round(mean(mu_ZVCV), 4), mse[2], round(mean(mu_CF), 4), mse[3]))

###########################################################
### Case 4: f(x) = exp(x) p(x) = 0.5 N(-2,1) + 0.5 N(2,1) ####
###########################################################

p <- 0.5
mu1 <- -2
mu2 <- 2
sd1 <- 1
sd2 <- 1

for (i in 1:rep){
  if(i%%10 == 0) print(i)
  samples <- as.matrix(mh.mcmc(start = -3, p, mu1, mu2, sd1, sd2, N+500, 1)[501:(N+500)], nrow = N, ncol = 1)
  grad <- as.matrix(apply(samples, 1, u, p, mu1, mu2, sd1, sd2), nrow = n_samp, ncol = 1)
  integrands <- as.matrix(exp(samples), nrow = n_samp, ncol = 1)
  mu_AM[i] <- mean(integrands)
  mu_ZVCV[i] <- zvcv(integrand = integrands, samples = samples, derivatives = grad)$expectation
  
  samples <- unique(samples)
  n_samp <- nrow(samples)
  grad <- as.matrix(apply(samples, 1, u, p, mu1, mu2, sd1, sd2), nrow = n_samp, ncol = 1)
  integrands <- as.matrix(exp(samples), nrow = n_samp, ncol = 1)
  
  mu_sCF[i] <- CF(integrands = integrands, samples = samples, derivatives = grad, steinOrder = 1, kernel_function = "prodsim", sigma = c(.1,1))$expectation
  ind <- sample.int(n_samp, n_samp/2, replace = FALSE)
  mu_CF[i] <- CF(integrands = integrands, samples = samples, derivatives = grad, steinOrder = 1, kernel_function = "prodsim", sigma = c(0.1,1), est_inds = ind)$expectation
  mu_ZVCV[i] <- zvcv(integrand = integrands, samples = samples, derivatives = grad)$expectation
}

pdf(file = "gaussian_mixture_mix2_exp.pdf", height = 5, width = 5)
boxplot(data.frame("AM" = mu_AM, "ZV-CV" = mu_ZVCV, "simplified CF" = mu_sCF), ylab = expression(mu(f)))
dev.off()

truth <- (exp(-1.5) + exp(2.5))/2
mse <- round(colMeans((cbind(mu_AM, mu_ZVCV, mu_sCF) - truth)^2),4)
print(paste(round(mean(mu_AM),4), mse[1], round(mean(mu_ZVCV), 4), mse[2], round(mean(mu_CF), 4), mse[3]))
