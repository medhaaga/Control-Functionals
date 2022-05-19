## Example
library(ZVCV)

n <- 200 # number of samples
d <- 10 # dimension of the state space
k <- 1 # dimension of range of f
rep <- 10
u <- function(x) -x # score vector (d/dx) log p(x)
f <- function(x){
  d <- length(x)
  return(1 + sin(pi*sum(x)/d))
} # integrand of interest
split = 0.5
m <- ceiling(n*split)
alpha = c(0.1, 1) # hyper-parameters
mu = 1 # true integral

mu_MC_rep <- rep(0, rep)
mu_CF_rep <- rep(0, rep)
mu_CFs_rep <- rep(0, rep)
mu_ZVCV_rep <- rep(0, rep)
mu_ZVCVs_rep <- rep(0, rep)

for (j in 1:rep){
  set.seed(j)
  samples <- matrix(rnorm(n*d), nrow = n, ncol = d) # random samples
  est_inds <- sample.int(n = n, size = m, replace = F)
  integrands <- matrix(0, n, k)
  derivatives <- matrix(0, n, d)
  for (i in 1:n){
    integrands[i,] <- f(samples[i,])
    derivatives[i,] <- u(samples[i,])
  } 
  
  mu_MC_rep[j] = colMeans(integrands) # classical Monte Carlo estimate
  CF_est = CF(integrand=integrands, samples=samples, derivatives=derivatives, steinOrder = 1, 
              est_inds = est_inds, kernel_function = "product", sigma = c(0.1, 1), )
  CFs_est = CF(integrand=integrands, samples=samples, derivatives=derivatives, steinOrder = 1, 
              kernel_function = "product", sigma = c(0.1, 1), )
 
  
  ZVCV_est = zvcv(integrand=integrands, samples=samples, derivatives=derivatives, est_inds = est_inds)
  ZVCVs_est = zvcv(integrand=integrands, samples=samples, derivatives=derivatives)
  
  mu_CF_rep[j] <- CF_est$expectation
  mu_CFs_rep[j] <- CFs_est$expectation
  mu_ZVCV_rep[j] <- ZVCV_est$expectation
  mu_ZVCVs_rep[j] <- ZVCVs_est$expectation
  
  if(j %% 10 == 0) print(j)
}

df <- data.frame("Monte Carlo" = mu_MC_rep, "CF" = mu_CF_rep, "CF-simplified" = mu_CFs_rep, 
                 "ZV-CV" = mu_ZVCV_rep, "ZVCV-simplified" = mu_ZVCVs_rep)
boxplot(df)
abline(h = 1)
