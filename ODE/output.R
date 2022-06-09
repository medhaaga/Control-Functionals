set.seed(1)
library(ZVCV)
source("functions.R")

y <- read.csv("y.csv", header = FALSE) ## observed data
colnames(y) <- NULL
rownames(y) <- NULL
y <- as.numeric(y)
n <- length(y) ## n = 11


m <- 30 ## no. of temperature intervals
temp <- (seq(0, m, 1)/m)^5  ## temperature intervals
N = 100 ## no. of MCMC samples

load(file = "Objects/100_sims.Rdata")
rep <- dim(theta_samp)[3]

tmp_pts <- seq(1,5,1)*20
mu_AM <- array(0, dim = c(length(temp), length(tmp_pts), rep))#matrix(0, nrow = length(temp), ncol = rep)
mu_CF <- array(0, dim = c(length(temp), length(tmp_pts), rep))#matrix(0, nrow = length(temp), ncol = rep)
mu_ZVCV <- array(0, dim = c(length(temp), length(tmp_pts), rep))
var_AM <- array(0, dim = c(length(temp), length(tmp_pts), rep))#matrix(0, nrow = length(temp), ncol = rep)
var_CF <- array(0, dim = c(length(temp), length(tmp_pts), rep))#matrix(0, nrow = length(temp), ncol = rep)
var_ZVCV <- array(0, dim = c(length(temp), length(tmp_pts), rep))

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

save(mu_AM, var_AM, mu_CF, var_CF, mu_ZVCV, var_ZVCV, file = "Objects/means_variances.Rdata")

load(file = "Objects/means_variances.Rdata")
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
marg_AM <- as.data.frame(marg_AM)
names(marg_AM) <-  c('20', '40', '60', '80', '100')
boxplot(marg_AM, ylim = range(marg_AM[,1], marg_CF[,1], marg_ZVCV[,1]), main = "Standard TI", xlab = "Num samples n", ylab = "Estimator Distribution", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
marg_ZVCV <- as.data.frame(marg_ZVCV)
names(marg_ZVCV) <-  c('20', '40', '60', '80', '100')
boxplot(marg_ZVCV, ylim = range(marg_AM[,1], marg_CF[,1], marg_ZVCV[,1]), main = "Controlled TI", xlab = "Num samples n", ylab = "Estimator Distribution", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
marg_CF <- as.data.frame(marg_CF)
names(marg_CF) <-  c('20', '40', '60', '80', '100')
boxplot(marg_CF, ylim = range(marg_AM[,1], marg_CF[,1], marg_ZVCV[,1]), main = "TI + Control Functionals", xlab = "Num samples n", ylab = "Estimator Distribution", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
dev.off()
