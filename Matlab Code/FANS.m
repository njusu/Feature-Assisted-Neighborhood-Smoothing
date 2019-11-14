function [P_hat] = FANS(A,X,lambda)
	% INPUT:
	% A: adjacency matrix
    % X: feature matrix
    % lambda: weight of distance based on feature, or chosen by CV
	% OUTPUT:
	% P_hat: estimated probaility matrix
	% COPYRIGHT:
	% Network estimation via graphon with nodefeatures
	% Yi Su, Raymond K.W. Wong, Thomas C.M. Lee
	
	N = size(A, 1);
    if N ~= size(X, 1)
        fprintf('WARNING!'); return;
    end
	h = sqrt(log(N)/N);
	
	% Compute disimilarity measures
	D0 = zeros(N, N);
	A_sq = A*A;
	for i = 1:(N-1)
		for j = (i+1):N
			temp = abs(A_sq(i, :) - A_sq(j, :));
            temp([i,j])=[];
            D0(i, j) = ( max(temp) +rand(1,1)/N )/N;  
			D0(j, i) = D0(i, j);
		end
    end
    
    s = pdist2(X,X);
    D = D0/max(D0(:)) + lambda * s/max(s(:)); 
    
    Kernel_mat = zeros(N,N);
	for i = 1:N
        Di = D(i,:); Di(i)=[];
		Kernel_mat(i, :) =  (D(i, :) <= quantile(Di, h));
        if(sum(Kernel_mat(i, :)) > 1)
            Kernel_mat(i, i) = 0;
        end
	end
	
	Kernel_mat = Kernel_mat ./ (sum(Kernel_mat, 2)*ones(1, N));  % Each row has been normalized(under L1)
	
	W_hat = Kernel_mat * A;
	P_hat = (W_hat + W_hat') / 2;
	
end