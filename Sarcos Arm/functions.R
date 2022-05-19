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

K_robot <- function(dist, theta){
  K <- matrix(0, nrow = nrow(dist), ncol = ncol(dist))
  for (i in 1:nrow(dist)){
    for (j in 1:ncol(dist)){
      K[i,j] <- theta[1]*exp(-dist[i,j]/(2*theta[2]^2))
    }
  }
  return(K)
}

