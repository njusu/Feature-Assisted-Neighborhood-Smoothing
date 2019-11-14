function y = g3(u,v,n)
y = 1-1/(1+exp(15*(0.8*abs(u-v))^(4/5)-0.1));