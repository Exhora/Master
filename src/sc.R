#clusterizacao espectral com medoids (PAM)
specClust <- function(W, k) {
  n <- ncol(W)
  ## Diagonal matrix with the degrees in the main diagonal
  D <- matrix(0, n, n)
  for (i in 1:n) {
    D[i,i] <- sum(abs(W[i,]))
  }
  ## Unnormalized Laplacian
  L <- D - W
  U <- (eigen(L)$vectors)[,((n-k+1):n)]
  C <- pam(x=U, k=k)
  return(C$clustering)
}
