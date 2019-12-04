# correlation between features and graphon slices
# Use Kendall tau correlation
library(pcaPP)
get_cor = function(d0, s){
  n = nrow(d0)
  #COR = rep(NA,n) # correlation for each individual
  #for(i in 1:n){ COR[i] = (cor(s[i,-i],d0[i,-i]))}
  diag(d0) = diag(s) = NA
  s_vec = as.vector(s); s_vec = s_vec[!is.na(s_vec)]
  d0_vec = as.vector(d0); d0_vec = d0_vec[!is.na(d0_vec)]
  #cor1 = cor(s_vec, d0_vec)
  #cor2 = cor.test(s_vec, d0_vec, method="spearman",exact = FALSE) # correlation for all entries
  kendalltau_cor = cor.fk(s_vec, d0_vec) # a fast estimation of kendall tau
  return(kendalltau_cor)
}
#get_cor(d0, s)

