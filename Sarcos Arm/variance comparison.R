load(file = "objects/means_n100.Rdata")


Ntest <- dim(mu_bar)[1]
sd_AM <- matrix(0, nrow = Ntest, ncol = 5)
sd_ZV <- matrix(0, nrow = Ntest, ncol = 5)
sd_CF <- matrix(0, nrow = Ntest, ncol = 5)
sd_CF_simp <- matrix(0, nrow = Ntest, ncol = 5)

for (i in 1:3){
  sd_AM[,i] <- apply(mu_bar[,,i], MARGIN = 1, FUN = sd)
  sd_ZV[,i] <- apply(mu_ZV[,,i], MARGIN = 1, FUN = sd)
  sd_CF[,i] <- apply(mu_CF[,,i], MARGIN = 1, FUN = sd)
  sd_CF_simp[,i] <- apply(mu_CF_simp[,,i], MARGIN = 1, FUN = sd)
}

load(file = "objects/means_n200.Rdata")

for (i in 1:2){
  sd_AM[,i+3] <- apply(mu_bar[,,i], MARGIN = 1, FUN = sd)
  sd_ZV[,i+3] <- apply(mu_ZV[,,i], MARGIN = 1, FUN = sd)
  sd_CF[,i+3] <- apply(mu_CF[,,i], MARGIN = 1, FUN = sd)
  sd_CF_simp[,i+3] <- apply(mu_CF_simp[,,i], MARGIN = 1, FUN = sd)
}


par(mfrow = c(1,3))
plot(sd_CF[,1], sd_AM[,1], xlim = c(0, 0.08), ylim = c(0, 0.08))
plot(sd_CF[,3], sd_AM[,3], xlim = c(0, 0.08), ylim = c(0, 0.08))
plot(sd_CF[,5], sd_AM[,5], xlim = c(0, 0.08), ylim = c(0, 0.08))

plot(sd_CF_simp[,1], sd_ZV[,1], xlim = c(0, 0.03), ylim = c(0, 0.03))
plot(sd_CF_simp[,3], sd_ZV[,3], xlim = c(0, 0.03), ylim = c(0, 0.03))
plot(sd_CF_simp[,5], sd_ZV[,5], xlim = c(0, 0.03), ylim = c(0, 0.03))

