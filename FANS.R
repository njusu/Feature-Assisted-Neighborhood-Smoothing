# FANS algorithm

FANS = function(A, X, lambda){
  # Estimate P_hat from A and X
  # Return P_hat
  
  get_d0 = function(A){
    n = nrow(A)
    d0 = matrix(0, n, n)
    A_sq = A%*%A
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        d0[i,j] = max( abs(A_sq[i,]-A_sq[j,])[-c(i,j)] +runif(1)/n )/n  #######
        d0[j,i] = d0[i,j]
      }
    }
    return(d0)
  }
  
  #s = as.matrix(dist(X, diag = TRUE)) # Euclidean dist
  get_s = function(X){
    n = nrow(X)
    s = matrix(0, n, n) # no sqrt
    X_sq = X%*%t(X)
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        s[i,j] = max( abs(X_sq[i,]-X_sq[j,])[-c(i,j)] +runif(1)/n )/ncol(X)  #######
        s[j,i] = s[i,j]
      }
    }
    return(s)
  }
  
  d0 = get_d0(A)
  s = get_s(as.matrix(X))
  
  d = d0/max(d0) + lambda*s/max(s)
  # Define neighborhood based on d(i,j)
  NB = matrix(NA,n,n)
  h = sqrt(log(n)/n)
  for(i in 1:n){
    NB[i,] = as.numeric(d[i,]<=quantile(d[i,-i],h))  ############## < ###############
    if(sum(NB[i,])>1){NB[i,i]=0 }
  }
  #which(NB[1,]!=0) # Neighbor of node i
  
  # Get P_hat
  kernel = NB # Just use (unweighted) average
  # can use other kernel matrix, e.g., kernel = NB*s, where s is the disimilarity matrix based on features
  W = kernel%*%A / (diag(rowSums(kernel)) %*% matrix(1,n,n))
  #return(W)  # asymmetric
  return( (W+t(W))/2 )  # symmetric
}
