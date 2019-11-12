# Cross-validation for FANS

CV = function(A, X, LAMBDA, ITER){
  #set.seed(Sys.time())
  #set.seed(1)
  nlambda = length(LAMBDA)
  PE = matrix(0,ITER,nlambda)  #PE is mean |P_hat-A|
  
  n = dim(A)[1]
  # Calculate distance
  d0 = matrix(0, n, n)
  A_sq = A%*%A
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      d0[i,j] = sqrt( max( abs(A_sq[i,]-A_sq[j,])[-c(i,j)] +runif(1)/n )/n )
      d0[j,i] = d0[i,j]
    }
  }
  s = as.matrix(dist(X, diag = TRUE))
  
  for(iter in 1:ITER){
    # Split data
    ind_valid = sort(sample(n,n%/%10)) # can sort or not sort, order doesn't matter
    ind_train = (1:n)[-ind_valid]
    A_train = A[ind_train,ind_train]
    X_train = as.matrix(X[ind_train,])
    X_valid = as.matrix(X[ind_valid,])
    
    # Fit training model
    n_train = length(ind_train)
    d0_train = matrix(0, n_train, n_train)
    A_sq_train = A_train%*%A_train
    for(i in 1:(n_train-1)){
      for(j in (i+1):n_train){
        d0_train[i,j] = sqrt( max( abs(A_sq_train[i,]-A_sq_train[j,])[-c(i,j)] +runif(1)/n )/n_train )   ##################
        d0_train[j,i] = d0_train[i,j]
      }
    }
    
    s_train = s[ind_train,ind_train]
    
    for(k in 1:nlambda){
      lambda = LAMBDA[k]
      d_train = d0_train/max(d0_train) + lambda*s_train/max(s_train)
      # Define neighborhood based on d(i,j)
      NB_train = matrix(NA,n_train,n_train)
      h_train = 1 * sqrt(log(n_train)/n_train)  # C=1
      for(i in 1:n_train){
        NB_train[i,] = as.numeric(d_train[i,]<=quantile(d_train[i,-i],h_train)) 
        if(sum(NB_train[i,])>1){NB_train[i,i]=0 }
      }
      #which(NB_train[1,]!=0) # Neighbor of node i
      
      # Get P_hat
      kernel_train = NB_train # Just use (unweighted) average
      W_train = kernel_train%*%A_train / (diag(rowSums(kernel_train)) %*% matrix(1,n_train,n_train))
      P_hat_train = (W_train+t(W_train))/2
      
      get_nearest = function(i){ which.min(s[i,ind_train]) } # i is from testing set, not in training set
      # Return an index in the training set, from 1 to 450.
      index = sapply(ind_valid, get_nearest)
      PE[iter,k] = mean(abs(P_hat_train[index,]-A[ind_valid,ind_train])) # abs error
    }
    
    if(iter%%1==0) print(iter)
  }
  return( colMeans(PE) )
}

#LAMBDA = c(0, 0.05*(1:60))
#res = CV(A,X,LAMBDA,20)
#opt = which.min(res)
#lambda_opt = LAMBDA[opt]