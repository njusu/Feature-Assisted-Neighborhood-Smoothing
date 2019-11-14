function y = g1(u,v,n)
K = floor(log(n));
for k=1:K
  if u>(k-1)/K && u<k/K && v>(k-1)/K && v<k/K
      y = k/(K+1); return;
  end
end
y = 0.3/(K+1);
end
