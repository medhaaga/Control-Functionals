set.seed(1)
library(ZVCV)
source("functions.R")

y <- read.csv("y.csv", header = FALSE) ## observed data
colnames(y) <- NULL
rownames(y) <- NULL
y <- as.numeric(y)
n <- length(y) ## n = 11

m <- 30 ## no. of temperatiure intervals
temp <- (seq(0, m, 1)/m)^5  ## temperatiure intervals
N = 100## no. of MCMC samples
rep <- 100

theta_samp <- array(0, dim = c(length(temp), N, rep))
f_theta <- array(0, dim = c(length(temp), N, rep)) ## f(theta) = log likelihood: log(p(y|theta))
u <- array(0, dim = c(length(temp), N, rep))       ## gradients of power density

## don't run this loop, objects are saved
for (itr in 1:rep){
  for (i in 1:length(temp)){
    
    print(paste("Replication = ", itr, "Temperature iter = ", i))
    theta_samp[i,,itr] <- MH(t = temp[i], N, sigma = 0.1, y)
    
    for (j in 1:N){
      x_a <- RK(theta_samp[i,j,itr])
      x_a_prime <- RK_deriv(theta_samp[i,j,itr], k=0.01)
      f_theta[i,j,itr] = (loglikelihood(y, theta_samp[i,j,itr], sigma = 0.1, x_a))
      u[i,j,itr] <- power_grad(theta_samp[i,j,itr], k = 0.001, t = temp[i], y, sigma = 0.1, x_a)
    }
  }
}

save(theta_samp, f_theta, u, file = "Objects/100_sims.Rdata")
