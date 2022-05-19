set.seed(1)
library(ZVCV)
source("functions.R")
set.seed(1)

Xtrain <- read.csv("Xtrain.csv", header = F)
Xtest <- read.csv("Xtest.csv", header = F)
Ytrain <- read.csv("Ytrain.csv", header = F)
Ytest <- read.csv("Ytest.csv", header = F)

muX <- colMeans(Xtrain)
sigmaX <- apply(Xtrain, 2, sd)
muY <- colMeans(Ytrain)
sigmaY <- apply(Ytrain, 2, sd)

Xtrain <- as.data.frame(scale(Xtrain))
Ytrain <- as.data.frame(scale(Ytrain))
Xtest <- as.data.frame(scale(Xtest, center = muX, scale = sigmaX))
Ytest <- as.data.frame(scale(Ytest, center = muY, scale = sigmaY))

Ntrain <- nrow(Xtrain)
N <- 1000
Nprime <- 100
ix <- sample.int(n = Ntrain, size = N, replace = F)
Xtrain <- Xtrain[ix,]
Ytrain <- as.matrix(Ytrain[ix,], ncol = 1)


ix <- sample.int(n = N, size = Nprime, replace = F)
Xprime = Xtrain[ix,]
Yprime = Ytrain[ix,]
Ntest <- nrow(Xtest)

split = 0.5

load(file = "distance_matrices.Rdata")
ns <- c(150, 200)
its <- 10
mu_bar <- array(0, dim = c(Ntest, its, length(ns)))
mu_CF <- array(0, dim = c(Ntest, its, length(ns)))
mu_CF_simp <- array(0, dim = c(Ntest, its, length(ns)))
mu_ZV <- array(0, dim = c(Ntest, its, length(ns)))

start <- Sys.time()

for (k in 1:length(ns)){
  print(paste("Running for N =", ns[k]))
  
  # Underlying distribution
  a <- 25
  b <- 0.04
  c <- 25
  d <- 0.04 # max cov / length scale
  u <- function(theta, a, b, c, d) return(c((a-1)/theta[1] - 1/b, (c-1)/theta[2] - 1/d)) # score
  n <- ns[k] # number of samples from pi
  m <- ceiling(n*split)
  est_inds <- sample.int(n = n, size = m, replace = F)
  
  theta = array(0, dim = c(n, 2, its))
  for (j in 1:its){
    theta[,1,j] <- rgamma(n, shape = a, scale = b)
    theta[,2,j] <- rgamma(n, shape = c, scale = d)
  }
  derivatives = array(0, dim = c(n, 2, its))
  for (j in 1:its){
    derivatives[,,j] <- t(apply(theta[,,j], MARGIN = 1, FUN = u, a, b, c, d))
  }
  
  
  # Target function 
  print("Evaluating functional values")
  f_vals = array(0, dim = c(Ntest,n,its))
  #h = waitbar(0,'evaluating function values');
  for (j in 1:its){
    for (i in 1:n){
      print(paste("Functional evaluation rep = ", j, "; iter = ", i))
      K_star_Nprime <- K_robot(D_star_Nprime, theta[i,,j])
      K_Nprime <- K_robot(D_Nprime, theta[i,,j])
      K_Nprime_N <- K_robot(D_Nprime_N, theta[i,,j])
      K_N_Nprime <- t(K_Nprime_N)
      sigma <- 0.1
      f_vals[,i,j] = K_star_Nprime %*% solve(K_Nprime_N %*% K_N_Nprime + sigma^2*K_Nprime) %*% (K_Nprime_N %*% Ytrain)
    }
  }
  save(f_vals, file = paste("objects/integrands_n", n, ".Rdata", sep = ""))
  
  # Estimation
  l <- 1
  sigma <- 0
  
  for (j in 1:its){
    for (i in 1:Ntest){
      mu_bar[i,j,k] <- mean(f_vals[i,,j])
      mu_ZV[i,j,k] <- zvcv(integrand = as.matrix(f_vals[i,,j]), samples = theta[,,j], derivatives = derivatives[,,j])$expectation
      mu_CF[i,j,k] <- CF(integrands = as.matrix(f_vals[i,,j]), samples = theta[,,j], derivatives = derivatives[,,j], steinOrder = 1, 
                         est_inds = est_inds, kernel_function = "product", sigma = c(0.1, 1))$expectation
      mu_CF_simp[i,j,k] <- CF(integrands = as.matrix(f_vals[i,,j]), samples = theta[,,j], derivatives = derivatives[,,j], steinOrder = 1, 
                              kernel_function = "product", sigma = c(0.1, 1))$expectation
    }
  }
  
  save(mu_bar, mu_ZV, mu_CF, mu_CF_simp, file = "objects/means_n200.Rdata")
}
stop <- Sys.time()
print(stop - start)


