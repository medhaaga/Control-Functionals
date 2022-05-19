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


m <- 30 ## no. of temperatiure intervals
temp <- (seq(0, m, 1)/m)^5  ## temperatiure intervals
N = 100## no. of MCMC samples
rep <- 10

theta_samp <- array(0, dim = c(length(temp), N, rep))
f_theta <- array(0, dim = c(length(temp), N, rep)) ## f(theta) = log likelihood: log(p(y|theta))
u <- array(0, dim = c(length(temp), N, rep))       ## gradients of power density

tmp_pts <- seq(1,5,1)*20
mu_AM <- array(0, dim = c(length(temp), length(tmp_pts), rep))#matrix(0, nrow = length(temp), ncol = rep)
mu_CF <- array(0, dim = c(length(temp), length(tmp_pts), rep))#matrix(0, nrow = length(temp), ncol = rep)
mu_ZVCV <- array(0, dim = c(length(temp), length(tmp_pts), rep))
var_AM <- array(0, dim = c(length(temp), length(tmp_pts), rep))#matrix(0, nrow = length(temp), ncol = rep)
var_CF <- array(0, dim = c(length(temp), length(tmp_pts), rep))#matrix(0, nrow = length(temp), ncol = rep)
var_ZVCV <- array(0, dim = c(length(temp), length(tmp_pts), rep))

## don't run this loop, objects are saved
for (itr in 1:rep){
  for (i in 1:length(temp)){
    
    print(paste("Replication = ", itr, "Temperature iter = ", i))
    theta_samp[i,,itr] <- MH(t = temp[i], N, sigma = 0.1, y)
    
    for (j in 1:N){
      x_a <- RK(theta_samp[i,j,itr])
      x_a_prime <- RK_deriv(theta_samp[i,j,itr], k=0.01)
      f_theta[i,j,itr] = (loglikelihood(y, theta_samp[i,j,itr], sigma = 0.1, x_a))
      u[i,j,itr] <- power_grad(theta_samp[i,j,itr], k = 0.01, t = temp[i], y, sigma = 0.1, x_a)
    }
  }
}

save(theta_samp, f_theta, u, file = "Objects/100_sims.Rdata")

for (itr in 1:rep){
  print(paste("Replication = ", itr))
  for (i in 1:length(temp)){
    for (tmp in 1:length(tmp_pts)){
      mu_AM[i,tmp,itr] <- mean(f_theta[i,(1:tmp_pts[tmp]),itr])
      var_AM[i,tmp,itr] <- var(f_theta[i,(1:tmp_pts[tmp]),itr])
      #est_ind <- sample.int(tmp_pts[tmp], tmp_pts[tmp]/2)
      mu_CF[i,tmp,itr] <- CF(integrands = as.matrix(f_theta[i,(1:tmp_pts[tmp]),itr]), samples = theta_samp[i,(1:tmp_pts[tmp]),itr], derivatives = u[i,(1:tmp_pts[tmp]),itr], steinOrder = 1,
                             kernel_function = "product", sigma = c(0.1, 3))$expectation
      var_CF[i,tmp,itr] <- CF(integrands = as.matrix((f_theta[i,(1:tmp_pts[tmp]),itr] - mu_CF[i,tmp,itr])^2), samples = theta_samp[i,(1:tmp_pts[tmp]),itr], derivatives = u[i,(1:tmp_pts[tmp]),itr], steinOrder = 1,
                              kernel_function = "product", sigma = c(0.1, 3))$expectation
      
      mu_ZVCV[i,tmp,itr] <- zvcv(integrand = (f_theta[i,(1:tmp_pts[tmp]),itr]), samples = theta_samp[i,(1:tmp_pts[tmp]),itr], derivatives = u[i,(1:tmp_pts[tmp]),itr], options = list(nfolds = 5))$expectation
      var_ZVCV[i,tmp,itr] <- zvcv(integrand = ((f_theta[i,(1:tmp_pts[tmp]),itr] - mu_ZVCV[i,tmp,itr])^2), samples = theta_samp[i,(1:tmp_pts[tmp]),itr], derivatives = u[i,(1:tmp_pts[tmp]),itr], options = list(nfolds = 5))$expectation
    }
  }
}


marg_AM <- matrix(0, nrow = rep, ncol = length(tmp_pts))
marg_CF <- matrix(0, nrow = rep, ncol = length(tmp_pts))
marg_ZVCV <- matrix(0, nrow = rep, ncol = length(tmp_pts))



## Calculate p(y) using 4-quadrature
for (itr in 1:rep){
  for (tmp in 1:length(tmp_pts)){
    for (i in 1:m){
      temp_diff <- (temp[i+1] - temp[i])
      marg_AM[itr, tmp] <- marg_AM[itr, tmp] + (temp_diff*((mu_AM[i+1,tmp,itr] + mu_AM[i,tmp,itr])/2)) - (temp_diff^2)*((var_AM[i+1,tmp,itr] - var_AM[i,tmp,itr])/12)
      marg_CF[itr, tmp] <- marg_CF[itr, tmp] + (temp_diff*((mu_CF[i+1,tmp,itr] + mu_CF[i,tmp,itr])/2)) - (temp_diff^2)*((var_CF[i+1,tmp,itr] - var_CF[i,tmp,itr])/12)
      marg_ZVCV[itr, tmp] <- marg_ZVCV[itr, tmp] + (temp_diff*((mu_ZVCV[i+1,tmp,itr] + mu_ZVCV[i,tmp,itr])/2)) - (temp_diff^2)*((var_ZVCV[i+1,tmp,itr] - var_ZVCV[i,tmp,itr])/12)
    }
  }
}

pdf(file = "Figures/boxplots.pdf", width = 15, height = 5)
par(mfrow = c(1,3))
mar_AM <- as.data.frame(marg_AM)
names(mar_AM) <-  c('20', '40', '60', '80', '100')
boxplot(mar_AM, ylim = range(mar_AM[,1], mar_CF[,1], mar_ZVCV[,1]), main = "Standard TI", xlab = "Num samples n", ylab = "Estimator Distribution")
mar_ZVCV <- as.data.frame(marg_ZVCV)
names(mar_ZVCV) <-  c('20', '40', '60', '80', '100')
boxplot(mar_ZVCV, ylim = range(mar_AM[,1], mar_CF[,1], mar_ZVCV[,1]), main = "Controlled TI", xlab = "Num samples n", ylab = "Estimator Distribution")
mar_CF <- as.data.frame(marg_CF)
names(mar_CF) <-  c('20', '40', '60', '80', '100')
boxplot(mar_CF, ylim = range(mar_AM[,1], mar_CF[,1], mar_ZVCV[,1]), main = "TI + Control Functionals", xlab = "Num samples n", ylab = "Estimator Distribution")
dev.off()
