## Example
set.seed(1)
library(ZVCV)

u <- function(x) -x # score vector (d/dx) log p(x)
f <- function(x){
  d <- length(x)
  return(1 + sin(pi*sum(x)/d))
} # integrand of interest


n <- 200 # number of samples
k <- 1 # dimension of range of f
rep <- 100
split = 0.5
m <- ceiling(n*split)
alpha = c(0.1, 1) # hyper-parameters
mu = 1 # true integral

time_rep <- seq(2,10,1)*10
times <- length(time_rep)
mu_MC_rep <- matrix(0, nrow = rep, ncol = times)
mu_CF_rep <- matrix(0, nrow = rep, ncol = times)
mu_sCF_rep <- matrix(0, nrow = rep, ncol = times)
mu_ZVCV_rep <- matrix(0, nrow = rep, ncol = times)


###########################################
############ d=1 ##########################
###########################################

d <- 1 # dimension of the state space

for (j in 1:rep){
  set.seed(j)
  samples <- matrix(rnorm(n*d), nrow = n, ncol = d) # random samples
  integrands <- matrix(0, n, k)
  derivatives <- matrix(0, n, d)
  for (i in 1:n){
    integrands[i,] <- f(samples[i,])
    derivatives[i,] <- u(samples[i,])
  } 
  for (t in 1:(times)){
    n_samp <- time_rep[t]
    ind <- sample.int(n_samp, n_samp/2, replace=FALSE)
    mu_MC_rep[j,t] <- mean(integrands[1:n_samp,]) 
    mu_sCF_rep[j,t] <- CF(integrand=integrands[1:n_samp,], samples=samples[1:n_samp,], derivatives=derivatives[1:n_samp,], steinOrder = 1, 
                         kernel_function = "product", sigma = c(0.1, 1))$expectation
    mu_CF_rep[j,t] <- CF(integrand=integrands[1:n_samp,], samples=samples[1:n_samp,], derivatives=derivatives[1:n_samp,], steinOrder = 1, 
                         kernel_function = "product", sigma = c(0.1, 1), est_inds = ind)$expectation
    mu_ZVCV_rep[j,t] <- zvcv(integrand=integrands[1:n_samp,], samples=samples[1:n_samp,], derivatives=derivatives[1:n_samp,])$expectation
  }
  if(j %% 10 == 0) print(j)
}

colnames(mu_MC_rep) <- time_rep
colnames(mu_CF_rep) <- time_rep
colnames(mu_sCF_rep) <- time_rep
colnames(mu_ZVCV_rep) <- time_rep

save(mu_MC_rep, mu_CF_rep, mu_sCF_rep, mu_ZVCV_rep, file = "objects/d1.Rdata")

############################################################################

load(file = "objects/d1.Rdata")
pdf(file = "objects/gaussian_d1.pdf", width = 12, height = 10)
par(mfrow = c(2,2), mar=c(5,6,4,1)+.1)
boxplot(as.data.frame(mu_MC_rep), ylim = range(mu_MC_rep, mu_ZVCV_rep), main = "Arithmetic means", xlab = "Number of samples", ylab = expression(mu(f)), cex.lab = 1.5, cex.main = 1.5, cex.axis=1)
boxplot(as.data.frame(mu_ZVCV_rep), ylim = range(mu_MC_rep, mu_ZVCV_rep), main = "Control variates", xlab = "Number of samples", ylab = expression(mu(f)), cex.lab = 1.5, cex.main = 1.5, cex.axis=1)
boxplot(as.data.frame(mu_CF_rep), ylim = range(mu_MC_rep, mu_ZVCV_rep), main = "Control functionals", xlab = "Number of samples", ylab = expression(mu(f)), cex.lab = 1.5, cex.main = 1.5, cex.axis=1)
boxplot(as.data.frame(mu_sCF_rep), ylim = range(mu_MC_rep, mu_ZVCV_rep), main = "Simplified CF", xlab = "Number of samples", ylab = expression(mu(f)), cex.lab = 1.5, cex.main = 1.5, cex.axis=1)
dev.off()

var_AM <- time_rep*colMeans((mu_MC_rep - 1)^2)
var_ZVCV <- time_rep*colMeans((mu_ZVCV_rep - 1)^2)
var_CF <- time_rep*colMeans((mu_CF_rep - 1)^2)
var_sCF <- time_rep*colMeans((mu_sCF_rep - 1)^2)

pdf(file = "objects/MSE.pdf", width = 5, height = 5)
par(mfrow = c(1,1))
plot(time_rep, var_AM, type = "b", ylim = range(var_AM, var_ZVCV, var_CF, var_sCF), pch =19, xlab = "Number of samples n", ylab = "MSE x n")
lines(time_rep, var_ZVCV, type = "b", col = 2, lty=2, pch = 18, lwd=2)
lines(time_rep, var_CF, type = "b", col  =3, lty=3, pch =17, lwd=2)
lines(time_rep, var_sCF, type = "b", col = 4, lty=4, pch = 16, lwd=2)
legend("topright", legend = c("AM", "ZV-CV", "CF", "simplified CF"), col = 1:4, lty = 1:4, lwd = c(1,2,2,2), pch = seq(19,16), cex = .7)
dev.off()

###########################################
############ d=10 #########################
###########################################

d <- 10 # dimension of the state space

for (j in 1:rep){
  set.seed(j)
  samples <- matrix(rnorm(n*d), nrow = n, ncol = d) # random samples
  integrands <- matrix(0, n, k)
  derivatives <- matrix(0, n, d)
  for (i in 1:n){
    integrands[i,] <- f(samples[i,])
    derivatives[i,] <- u(samples[i,])
  } 
  for (t in 1:(times)){
    n_samp <- time_rep[t]
    ind <- sample.int(n_samp, n_samp/2, replace=FALSE)
    mu_MC_rep[j,t] <- mean(integrands[1:n_samp,]) 
    mu_sCF_rep[j,t] <- CF(integrand=integrands[1:n_samp,], samples=samples[1:n_samp,], derivatives=derivatives[1:n_samp,], steinOrder = 1, 
                          kernel_function = "product", sigma = c(0.1, 1))$expectation
    mu_CF_rep[j,t] <- CF(integrand=integrands[1:n_samp,], samples=samples[1:n_samp,], derivatives=derivatives[1:n_samp,], steinOrder = 1, 
                         kernel_function = "product", sigma = c(0.1, 1), est_inds = ind)$expectation
    mu_ZVCV_rep[j,t] <- zvcv(integrand=integrands[1:n_samp,], samples=samples[1:n_samp,], derivatives=derivatives[1:n_samp,])$expectation
  }
  if(j %% 10 == 0) print(j)
}

colnames(mu_MC_rep) <- time_rep
colnames(mu_CF_rep) <- time_rep
colnames(mu_sCF_rep) <- time_rep
colnames(mu_ZVCV_rep) <- time_rep

save(mu_MC_rep, mu_CF_rep, mu_sCF_rep, mu_ZVCV_rep, file = "objects/d10.Rdata")

############################################################################

load(file = "objects/d10.Rdata")
pdf(file = "objects/gaussian_d10.pdf", width = 12, height = 10)
par(mfrow = c(2,2), mar=c(5,6,4,1)+.1)
boxplot(as.data.frame(mu_MC_rep), ylim = range(mu_MC_rep, mu_ZVCV_rep), main = "Arithmetic means", xlab = "Number of samples", ylab = expression(mu(f)), cex.lab = 1.5, cex.main = 1.5, cex.axis=1)
boxplot(as.data.frame(mu_ZVCV_rep), ylim = range(mu_MC_rep, mu_ZVCV_rep), main = "Control variates", xlab = "Number of samples", ylab = expression(mu(f)), cex.lab = 1.5, cex.main = 1.5, cex.axis=1)
boxplot(as.data.frame(mu_CF_rep), ylim = range(mu_MC_rep, mu_ZVCV_rep), main = "Control functionals", xlab = "Number of samples", ylab = expression(mu(f)), cex.lab = 1.5, cex.main = 1.5, cex.axis=1)
boxplot(as.data.frame(mu_sCF_rep), ylim = range(mu_MC_rep, mu_ZVCV_rep), main = "Simplified CF", xlab = "Number of samples", ylab = expression(mu(f)), cex.lab = 1.5, cex.main = 1.5, cex.axis=1)
dev.off()


