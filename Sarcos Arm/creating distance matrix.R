distance <- function(X1, X2){
  D <- matrix(0, nrow = nrow(X1), ncol = nrow(X2))
  for (i in 1:nrow(X1)){
    for (j in 1:nrow(X2)){
      print(i)
      D[i,j] <- crossprod(as.numeric(X1[i,] - X2[j,]))
    }
  }
  return(D)
}
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


D_star_Nprime <- distance(Xtest, Xprime)
D_Nprime <- distance(Xprime, Xprime)
D_Nprime_N <- distance(Xprime, Xtrain)
save(D_star_Nprime, D_Nprime, D_Nprime_N, file = "objects/distance_matrices.Rdata")
