# Real data (Figures related to real data)
setwd("~/Google Drive/Research/Review IEEE")

library(fields)
source('FANS.R')
source('CV_for_FANS.R')

######### ADD HEALTH #########
ADD = as.matrix(read.table("~/Google Drive/Research/Final Results and Codes/Real data/comm10.dat"))
ADD_X = as.matrix(read.table("~/Google Drive/Research/Final Results and Codes/Real data/comm10_att.dat"))
n = nrow(ADD_X)
ADD_A = matrix(0,n,n)
for(i in 1:nrow(ADD)){
  ADD_A[ADD[i,1],ADD[i,2]]=ADD[i,3]
}
index_rm = which(rowSums(ADD_A)==0|ADD_X[,1]==0|ADD_X[,2]==0|ADD_X[,3]==0)
A = ADD_A[-index_rm,-index_rm]
A = 1*(A+t(A)>0)
X_raw = ADD_X[-index_rm,1:3]
n = nrow(X_raw)

# Convert feature 2 to an nx6 matrix
f2 = matrix(0, n, 6)
for(i in 1:n){
  f2[i, X_raw[i,2]+1] = 1
}
X = cbind(X_raw[,1],X_raw[,3],f2)

# sort by grade
dg_index = order(X[,2],decreasing = TRUE)
A = A[dg_index,dg_index]
X = X[dg_index,]
image(A, col=c('white','black'))

# Results
# CV
LAMBDA = c(0, 0.01*(1:50))
res = CV(A,X,LAMBDA,10)
opt = which.min(res)
lambda_opt = LAMBDA[opt]  # 0.12, chosen by CV

# Fit
t0 = proc.time()
fit = FANS(A,X,0.1)
proc.time() - t0
image(fit)
image.plot(fit, main = "FANS, lambda = 0.1", nlevel = 1000, zlim = c(0, 0.1))
# then sort by gender
dg_index = order(X[,1],decreasing = TRUE)
image.plot(fit[dg_index,dg_index], main = "FANS, lambda = 0.1, sorted by grade and gender", nlevel = 1000, zlim = c(0, 0.1))
proportion = sum(X[,1]==2)/nrow(A)
abline(v=proportion,col='yellow', lwd = 2, lty = 2)
abline(h=proportion,col='yellow', lwd = 2, lty = 2)

# Model Comparison
# AUC plot
library(AUC)
library(pROC)

# Write table as input to Matlab. Other methods are implemented in Matlab.
write.table(A, 'comm10_A.txt', row.names=FALSE, col.names=FALSE)

fit1 = FANS(A,X,0)
image.plot(fit1, main = "FANS, lambda = 0", nlevel = 1000, zlim = c(0, 0.1))
roc_obj_fans0 <- roc(as.factor(as.vector(A[row(A)!=col(A)])), as.vector(fit1[row(fit1)!=col(fit1)]))
auc(roc_obj_fans0) # 0.8867

fit2 = FANS(A,X,0.1)
image.plot(fit2, main = "FANS, lambda = 0.1", nlevel = 1000, zlim = c(0, 0.1))
roc_obj_fans01 <- roc(as.factor(as.vector(A[row(A)!=col(A)])), as.vector(fit2[row(fit1)!=col(fit2)]))
auc(roc_obj_fans01) # 0.9231

P_NBS = as.matrix(read.table('~/Google Drive/Research/Final Results and Codes/Real data/Model Comparison/comm10_NBS.txt',sep=','))
image.plot(P_NBS, main = "NBS", nlevel = 1000, zlim = c(0, 0.1))
roc_obj_nbs <- roc(as.factor(as.vector(A[row(A)!=col(A)])), as.vector(P_NBS[row(P_NBS)!=col(P_NBS)]))
auc(roc_obj_nbs) # 0.8247

P_SAS = as.matrix(read.table('~/Google Drive/Research/Final Results and Codes/Real data/Model Comparison/comm10_SAS.txt',sep=','))
# actually max(P_SAS)<0.1
image.plot(P_SAS, main = "SAS", nlevel = 1000, zlim = c(0, 0.1))
roc_obj_sas <- roc(as.factor(as.vector(A[row(A)!=col(A)])), as.vector(P_SAS[row(P_SAS)!=col(P_SAS)]))
auc(roc_obj_sas) # 0.7068

P_SAS_sorted = as.matrix(read.table('~/Google Drive/Research/Final Results and Codes/Real data/Model Comparison/comm10_SAS_sorted.txt',sep=','))
image.plot(P_SAS_sorted, main = "SAS (sorted by node degree)", nlevel = 1000, zlim = c(0, 0.1))

P_USVT = as.matrix(read.table('~/Google Drive/Research/Final Results and Codes/Real data/Model Comparison/comm10_USVT.txt',sep=','))
P_USVT[P_USVT>0.1] = 0.1
image.plot(P_USVT, main = "USVT", nlevel = 1000, zlim = c(0, 0.1))
roc_obj_usvt <- roc(as.factor(as.vector(A[row(A)!=col(A)])), as.vector(P_USVT[row(P_USVT)!=col(P_USVT)]))
auc(roc_obj_usvt) # 0.5

P_JCDC = as.matrix(read.table('~/Google Drive/Research/Final Results and Codes/Real data/Model Comparison/comm10_JCDC.txt'))
image.plot(P_JCDC, main = "JCDC, K = 6", nlevel = 1000, zlim = c(0, 0.1))
roc_obj_jcdc <- roc(as.factor(as.vector(A[row(A)!=col(A)])), as.vector(P_JCDC[row(P_JCDC)!=col(P_JCDC)]))
auc(roc_obj_jcdc) # 0.8776

# ROC curve comparison
plot(1-roc_obj_fans01$specificities,roc_obj_fans01$sensitivities, lwd = 2, col = 'blue',
     type='l',main = "ROC curve", xlab = "False Positive Rate", ylab = "True Positive Rate")
lines(1-roc_obj_fans0$specificities,roc_obj_fans0$sensitivities, lwd = 2, lty = 2, col = 'black')
lines(1-roc_obj_nbs$specificities,roc_obj_nbs$sensitivities, lwd = 2, lty = 3, col = 'red')
lines(1-roc_obj_sas$specificities,roc_obj_sas$sensitivities, lwd = 2, lty = 4, col = 'green')
lines(1-roc_obj_jcdc$specificities,roc_obj_jcdc$sensitivities, lwd = 2, lty = 5, col = 'orange')
lines(1-roc_obj_usvt$specificities,roc_obj_usvt$sensitivities, lwd = 2, lty = 6, col = 'pink')

# manually add border for legend
x_pos = 0.6; y_pos = 0.53
segments(x_pos, -0.5, x_pos, y_pos)
segments(x_pos, y_pos, 1.5, y_pos)
legend(x_pos, y_pos,
       legend = c('FANS (lambda = 0.1) \nAUC = 0.92',
                  'FANS (lambda = 0) \nAUC = 0.89',
                  'NBS \nAUC = 0.82',
                  'SAS \nAUC = 0.71',
                  'JCDC \nAUC = 0.88',
                  'USVT \nAUC = 0.5'),
       lwd = 2, lty = 1:6, col = c('blue', 'black', 'red', 'green', 'orange', 'pink'), cex = 1,
       x.intersp = 0.6, y.intersp = 0.5, bty = 'n')

