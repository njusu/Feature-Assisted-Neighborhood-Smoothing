# Study the effect of lambda, plot MSE against lambda (Figure 2, 4)
# Based on MSEplot_new.R in "Final Results and Codes" folder

# g0
g0 = function(u,v){ (u+v)/2}
# g1
g1 = function(u,v){
  K = as.integer(log(n))
  for(k in 1:K){
    if(u>(k-1)/K & u<k/K & v>(k-1)/K & v<k/K){ return(k/(K+1)) }
  }
  return(0.3/(K+1))
}
# g2
g2 = function(u,v){ sin(5*pi*(u+v-1)+1)/2+0.5 }
# g3
g3 = function(u,v){ 1-1/(1+exp(15*(0.8*abs(u-v))^(4/5)-0.1))}
# g4
g4 = function(u,v){ (u^2+v^2)/3*cos(1/(u^2+v^2))+0.15 }

# feature 0 (for illustrative example)
f0 = function(x){ cos(4*pi*u) } # periodic, 1,2,4
# feature 1
f1 = function(x){ 2*cos(2*pi*(1-x)^2)} # not monotone
# feature 2
f2 = function(x){ 10*x^2-12*x+5} # not monotone
# feature 3
f3 = function(x){ 2*cos(pi*x)} # monotone
# feature 4
f4 = function(x){ 2/3*qnorm(x)} # monotone


sigma = 0.3
n = 500
c = 1

LAMBDA = c(0, 0.05*(1:100))
maxRun = 20
MSE1 = MSE2 = MSE3 = MSE4 = matrix(NA, maxRun, length(LAMBDA))

g = g1
for(run in 1:maxRun){
  #set.seed(123) # If run!=1, delete this line
  u0 = v0 = runif(n)
  u = v = sort(u0)
  
  # Feature
  X = cbind(f1(u)/sd(f1(u))+rnorm(n, 0, sigma), f2(u)/sd(f2(u))+rnorm(n, 0, sigma),
            f3(u)/sd(f3(u))+rnorm(n, 0, sigma), f4(u)/sd(f4(u))+rnorm(n, 0, sigma))
  
  # True probability matrix P and observed adjacency matrix A
  P = A = matrix(0, n, n)
  for(i in 1:n){
    for(j in i:n){
      P[i,j] = g(u[i],u[j])
      A[i,j] = rbinom(1, 1, P[i,j])
      P[j,i] = P[i,j]
      A[j,i] = A[i,j]
    }
  }
  #diag(P)=diag(A)=0
  
  # Calculate distance (no sqrt)
  d0 = matrix(0, n, n)
  A_sq = A%*%A
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      d0[i,j] = max( abs(A_sq[i,]-A_sq[j,])[-c(i,j)] +runif(1)/n )/n  #######################
      d0[j,i] = d0[i,j]
    }
  }
  
  s = matrix(0, n, n) # no sqrt
  X_sq = X%*%t(X)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      s[i,j] = max( abs(X_sq[i,]-X_sq[j,])[-c(i,j)] +runif(1)/n )/ncol(X)
      s[j,i] = s[i,j]
    }
  }
  
  iter = 0
  for(lambda in LAMBDA){
    iter = iter+1
    # Define neighborhood based on d(i,j)
    d = d0/max(d0) + lambda*s/max(s)
    NB = matrix(NA,n,n)
    h = 1*sqrt(log(n)/n)
    
    for(i in 1:n){
      NB[i,] = as.numeric(d[i,]<=quantile(d[i,-i],h))
      if(sum(NB[i,])>1){NB[i,i]=0 }
    }

    # Get P_hat and MSE
    kernel = NB # Just use (unweighted) average
    # can use other kernel matrix, e.g., kernel = NB*sim1, where sim1 is the disimilarity matrix based on features
    W = kernel%*%A / (rowSums(kernel) %*% matrix(1,1,n))
    P_hat = (W+t(W))/2
    
    MSE1[run,iter] = mean((P_hat-P)^2)
    MSE2[run,iter] = mean(abs(P_hat-P))
    MSE3[run,iter] = mean((P_hat-A)^2)
    MSE4[run,iter] = mean(abs(P_hat-A))
  }
  print(run)
}

write.table(colMeans(MSE1),"g1_allf_noise0.3_MSE1.txt")
write.table(colMeans(MSE2),"g1_allf_noise0.3_MSE2.txt")
write.table(colMeans(MSE3),"g1_allf_noise0.3_MSE3.txt")
write.table(colMeans(MSE4),"g1_allf_noise0.3_MSE4.txt")


#####################
###### Results ######
#####################
plot(LAMBDA, colMeans(MSE1), main = "g1, n=200, allf with noise = 0 \n run=20, MSE vs lambda")

# Read server result
noise0 = read.table("g1_allf_noise0_MSE1.txt")$x
noise0.1 = read.table("g1_allf_noise0.1_MSE1.txt")$x
noise0.2 = read.table("g1_allf_noise0.2_MSE1.txt")$x
noise0.3 = read.table("g1_allf_noise0.3_MSE1.txt")$x
noise0.4 = read.table("g1_allf_noise0.4_MSE1.txt")$x
noise0.5 = read.table("g1_allf_noise0.5_MSE1.txt")$x

plot(LAMBDA, noise0, col = 1, type = 'l', lwd = 2, 
     main = 'MSE vs feature noise , g1 all features \n run=20 for each curve',
     xlab = 'lambda', ylab = 'MSE' )
points(LAMBDA, noise0.1, type = 'l', lwd = 2, col = 2)
points(LAMBDA, noise0.2, type = 'l', lwd = 2, col = 3)
points(LAMBDA, noise0.3, type = 'l', lwd = 2, col = 4)
points(LAMBDA, noise0.4, type = 'l', lwd = 2, col = 5)
points(LAMBDA, noise0.5, type = 'l', lwd = 2, col = 6)

legend('topright', legend = c(0,0.1,0.2,0.3,0.4,0.5),
       lwd = 2, col = 1:6, cex = 0.8, y.intersp = 0.7)

