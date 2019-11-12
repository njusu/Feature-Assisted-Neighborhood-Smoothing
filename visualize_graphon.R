# Visualize graphon functions and features (Figure 1, 2, 3, 11)

library('fields')

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


# feature 1
f1 = function(x){ 2*cos(2*pi*(1-x)^2)} # not monotone
# feature 2
f2 = function(x){ 10*x^2-12*x+5} # not monotone
# feature 3
f3 = function(x){ 2*cos(pi*x)} # monotone
# feature 4
f4 = function(x){ 2/3*qnorm(x)} # monotone

# Visualize features
x = seq(0.001, 0.999, 0.001)
f1x = sapply(x, f1)
f2x = sapply(x, f2)
f3x = sapply(x, f3)
f4x = sapply(x, f4)
plot(x, (f1x-mean(f1x))/sd(f1x), type='l', ylim = c(-3, 3), col=1, lwd=1.5,
     xlab='u', ylab='f(u)', main='Features used in simulation')
lines(x, (f2x-mean(f2x))/sd(f2x), col=2, lwd=1.5)
lines(x, (f3x-mean(f3x))/sd(f3x), col=3, lwd=1.5)
lines(x, (f4x-mean(f4x))/sd(f4x), col=4, lwd=1.5)
legend('bottomright', legend=c('feature 1', 'feature 2', 'feature 3', 'feature 4'),
       lty=1, lwd=1.5, col=1:4, y.intersp=0.6, x.intersp=0.5)

# Visualize graphon
n = 1000

g = g1
u0 = v0 = runif(n)
u = v = sort(u0)
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
image.plot(P, col=heat.colors(500), main="Graphon 1", zlim = c(0, 1), cex.main=1.5, cex.axis=1.5, legend.cex=1.5)
